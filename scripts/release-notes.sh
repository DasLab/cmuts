#!/usr/bin/env bash
#
# release-notes.sh -- print the changelog entry for one version.
#
# usage: scripts/release-notes.sh VERSION
#
# The entry is what CHANGELOG.md holds under `## [VERSION]`, up to the heading
# of the entry that follows it. The release workflow attaches this to the draft
# it creates, so that the release carries the same text as the changelog.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -euo pipefail

PROGRAM="release-notes.sh"

die() {
    printf '%s: %s\n' "$PROGRAM" "$1" >&2
    exit 1
}

[ $# -eq 1 ] || die "usage: $PROGRAM VERSION"

NOTES=$(awk -v version="$1" '
    $0 == "## [" version "]" || index($0, "## [" version "] ") == 1 { inside = 1; next }
    inside && /^## / { exit }
    inside { print }
' CHANGELOG.md)

# The blank line under the heading. The trailing ones go with the substitution
# above, which drops every newline at the end.
NOTES=$(printf '%s\n' "$NOTES" | sed '/./,$!d')

[ -n "$NOTES" ] || die "CHANGELOG.md holds no entry for $1"

printf '%s\n' "$NOTES"
