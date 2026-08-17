#!/usr/bin/env bash
# Best-effort, guarded database setup for Yggdrasil.
# Usage: bash setup_databases.sh <db_dir>
# Each tool's DB is large and tool-specific; this script creates the layout,
# runs whichever downloader is available on PATH, and otherwise prints the
# exact command to run inside that tool's conda env (created by the workflow).
set -euo pipefail
DB="${1:-databases}"
mkdir -p "$DB"/{checkv,vcontact3,phrog,pharokka,dram,iphop,virclust,viridic,phold}

have() { command -v "$1" >/dev/null 2>&1; }
note() { echo "[setup-databases] $*"; }

# CheckV
if [[ -f "$DB/checkv/.done" ]]; then note "CheckV DB present"; else
  if have checkv; then checkv download_database "$DB/checkv" && touch "$DB/checkv/.done";
  else note "CheckV: run inside its env ->  conda run -n <checkv-env> checkv download_database $DB/checkv"; fi
fi

# Pharokka / PHROG (pharokka bundles PHROG + mmseqs DBs)
if [[ -f "$DB/pharokka/.done" ]]; then note "Pharokka DB present"; else
  if have pharokka.py; then pharokka.py --databases -o "$DB/pharokka" && touch "$DB/pharokka/.done";
  else note "Pharokka: run inside its env ->  pharokka.py --databases -o $DB/pharokka  (then point PHROG at $DB/phrog)"; fi
fi

# vContact3 (downloads on first run to a default path; set VCONTACT3_DB to reuse)
if [[ -f "$DB/vcontact3/.done" ]]; then note "vContact3 DB present"; else
  note "vContact3: it fetches its DB on first run; set VCONTACT3_DB=$DB/vcontact3 or download from the vcontact3 data release into $DB/vcontact3, then touch $DB/vcontact3/.done";
fi

# iPHoP (optional). Full DB is ~300 GB (default `iphop download`); a ~9 GB test DB exists for wiring checks.
if [[ -f "$DB/iphop/.done" ]]; then note "iPHoP DB present"; else
  if have iphop; then iphop download --db_dir "$DB/iphop" --no_prompt && touch "$DB/iphop/.done";
  else note "iPHoP (optional, in-env): iphop download --db_dir $DB/iphop --no_prompt   # ~300GB; test wiring: add -dbv iPHoP_db_rw_1.4_for-test (~9GB)"; fi
fi

# DRAM-v (optional) - large DB; AMG detection here defaults to phold instead (below)
note "DRAM-v (optional, NOT used by default): DRAM-setup.py prepare_databases --output_dir $DB/dram  (large)"

# phold (optional; default AMG path) - structure-informed phage annotation DB (~16 GB).
# The amg_phold_db rule auto-installs it; this is here for manual setup / documentation.
if [[ -f "$DB/phold/.done" ]]; then note "phold DB present"; else
  note "phold (AMG, in-env): phold install -d $DB/phold --foldseek_gpu -t 8   # ~16GB; drop --foldseek_gpu on CPU-only"; fi

# VirClust (optional): clone standalone R scripts. Heavy R env (envs/virclust.yaml); run only when virclust.enabled.
if [[ -f "$DB/virclust/.done" ]]; then note "VirClust scripts present"; else
  if have git; then git clone --depth 1 https://github.com/CristinaMoraru/VirClust.git "$DB/virclust/repo" && touch "$DB/virclust/.done";
  else note "VirClust: git not found ->  git clone --depth 1 https://github.com/CristinaMoraru/VirClust.git $DB/virclust/repo"; fi
fi

# VIRIDIC (optional): clone standalone R scripts. Heavy R env (envs/viridic.yaml); run only when viridic.enabled.
if [[ -f "$DB/viridic/.done" ]]; then note "VIRIDIC scripts present"; else
  if have git; then git clone --depth 1 https://github.com/CristinaMoraru/VIRIDIC.git "$DB/viridic/repo" && touch "$DB/viridic/.done";
  else note "VIRIDIC: git not found ->  git clone --depth 1 https://github.com/CristinaMoraru/VIRIDIC.git $DB/viridic/repo"; fi
fi

note "Database layout ready under $DB (some steps above may be manual)."
