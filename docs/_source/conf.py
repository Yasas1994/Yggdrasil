# Sphinx configuration for Yggdrasil docs.
project = "yggdrasil"
copyright = "2026, Yggdrasil contributors"
author = "Yggdrasil contributors"
html_title = "yggdrasil"
release = "0.1.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_design",
    "sphinx_copybutton",
]
source_suffix = [".rst", ".md"]
myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath"]
myst_heading_anchors = 3

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
pygments_style = "catppuccin-latte"
pygments_dark_style = "catppuccin-mocha"
html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "light_css_variables": {
        "color-brand-primary": "#1e66f5",
        "color-brand-content": "#1e66f5",
    },
    "dark_css_variables": {
        "color-brand-primary": "#8aadf4",
        "color-brand-content": "#8aadf4",
    },
}
html_static_path = ["_static"]
