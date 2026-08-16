# conf.py -- configuration for the documentation site.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

project = "cmuts"
html_title = "cmuts"
author = "Hamish M. Blair"

# The first year is the one in LICENSE, which covers the documentation as well.
# %Y is the year the site is built, so the footer does not go stale.
copyright = "2024-%Y, Hamish M. Blair"

extensions = [
    "myst_parser",
    "sphinx_inline_tabs",
]

# colon_fence is for the ::: fences the install tabs are written with.
# attrs_block is for the {.class} written above a block, which MyST attaches to
# that block, or to the section where the block is a heading.
myst_enable_extensions = ["colon_fence", "attrs_block"]

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Where the site is served from, which the canonical link of each page needs.
# The organization's root site carries a custom domain, so a project page of
# the organization is served under that rather than under github.io.
html_baseurl = "https://daslab.stanford.edu/cmuts/"

# Copied to the root of the site. GitHub Pages runs the output through Jekyll
# otherwise, which drops the _static and _sources directories sphinx writes.
# The file is empty; its presence is the whole of it.
html_extra_path = [".nojekyll"]

html_theme_options = {
    "source_repository": "https://github.com/DasLab/cmuts",
    "source_branch": "main",
    "source_directory": "docs/",
}
