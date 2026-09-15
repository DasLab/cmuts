#!/usr/bin/env bash
#
# check-tag.sh -- confirm a release tag names the version the programs report.
#
# usage: scripts/check-tag.sh TAG
#
# The release takes its name from the tag and every tarball in it takes its
# name from include/version.h, so a tag disagreeing with the header publishes
# binaries which report a version other than the one the release announces.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -euo pipefail

PROGRAM="check-tag.sh"

die() {
    printf '%s: %s\n' "$PROGRAM" "$1" >&2
    exit 1
}

[ $# -eq 1 ] || die "usage: $PROGRAM TAG"

VERSION=$(sed -n 's/.*CMUTS_VERSION "\(.*\)".*/\1/p' include/version.h)

[ -n "$VERSION" ] || die "include/version.h holds no version"
[ "$1" = "v$VERSION" ] || die "tag $1 names a version other than $VERSION, which include/version.h holds"

printf '%s: %s matches include/version.h\n' "$PROGRAM" "$1"
