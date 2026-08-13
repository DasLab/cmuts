# conf.py -- configuration for the documentation site.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

project = "cmuts"
html_title = "cmuts"

extensions = [
    "myst_parser",
    "sphinx_inline_tabs",
]

# For the ::: fences the install tabs are written with.
myst_enable_extensions = ["colon_fence"]

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "source_repository": "https://github.com/hmblair/cmuts",
    "source_branch": "main",
    "source_directory": "docs/",
}
