"""Sphinx configuration for the NumTracer documentation site.

The site bundles two language APIs under one theme:

* **C++** — Doxygen renders the API to HTML; ``doxysphinx`` wraps that HTML
  into this Sphinx tree (see ``build.sh``). Output lives in ``doxygen/`` and is
  generated, not committed.

The Mathematica codegen front-end is documented narratively under *Under the Hood*
(``internals/codegen.md``); there is intentionally no Python API (unlike the upstream
DiFfRG site this is modelled on), so ``sphinx-autoapi`` and the NumPy intersphinx
mapping are absent.
"""

from pathlib import Path

HERE = Path(__file__).parent.resolve()
NUMTRACER_DIR = HERE.parent  # the numtracer/ project root

# -- Project information ------------------------------------------------------

project = "NumTracer"
author = "Franz R. Sattler"
copyright = "2026, Franz R. Sattler"

release = "0.1.0"
version = release

# -- General configuration ----------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx_design",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"

# Do not let Sphinx read the bundled Doxygen theme, the Doxygen config, or build
# artefacts. The generated doxygen/ rst (from doxysphinx) MUST stay visible, so
# it is not excluded.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "README.md",
    "doxygen-awesome-css",
    "requirements.txt",
    # Never scan a virtualenv that happens to live in the source tree.
    ".venv*",
    "**/site-packages/**",
]

# -- MyST ---------------------------------------------------------------------

myst_enable_extensions = [
    "dollarmath",   # $...$ and $$...$$
    "amsmath",      # \begin{aligned} ... \end{aligned}
    "deflist",
    "colon_fence",
    "attrs_inline",
    "fieldlist",
    "html_image",
]
myst_heading_anchors = 3

# TeX macros for the physics notation used across the narrative pages.
mathjax3_config = {
    "tex": {
        "macros": {
            # Feynman slash: \slashed{p} (also accepts the unbraced \slashed p).
            "slashed": ["{#1\\!\\!\\!/}", 1],
        }
    }
}

# Pygments cannot strictly lex some of the C++/template signatures shown in the
# narrative pages; the highlighter retries in relaxed mode and renders fine, so
# silence the cosmetic failures.
suppress_warnings = ["misc.highlighting_failure"]

# -- HTML output --------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "NumTracer"
html_show_sourcelink = False  # drop the "Show Source" link

html_static_path = ["_static"]
# Make the embedded Doxygen (doxygen-awesome) pages follow the site's dark theme.
html_css_files = ["doxygen-darkmode.css"]

# The pydata theme builds the top navbar from the root document's top-level
# toctree entries — these become the Getting Started / Tutorials / Under the Hood
# / C++ API tabs.
html_theme_options = {
    "navbar_align": "left",
    # Auto-expand the left ("section navigation") sidebar two levels deep so the
    # API tree isn't a single collapsed entry.
    "show_nav_level": 2,
    "navigation_depth": 4,
    "collapse_navigation": False,
    "use_edit_page_button": False,
    # Single sidebar: fold the in-page "On this page" TOC away (no right column).
    "secondary_sidebar_items": [],
}
