#!/usr/bin/env bash
#
# check-tag.sh -- confirm a release tag names the version the programs report
# and the version the citation gives.
#
# usage: scripts/check-tag.sh TAG
#
# The release takes its name from the tag and every tarball in it takes its
# name from include/version.h, so a tag disagreeing with the header publishes
# binaries which report a version other than the one the release announces.
# CITATION.cff names the version and release date a citation of the software
# carries, and Zenodo reads it when it archives the tag, so a citation file
# left at an older version misnames the release everywhere it is cited.
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

CITED=$(sed -n 's/^version: *//p' CITATION.cff)

[ -n "$VERSION" ] || die "include/version.h holds no version"
[ "$1" = "v$VERSION" ] || die "tag $1 names a version other than $VERSION, which include/version.h holds"
[ -n "$CITED" ] || die "CITATION.cff holds no version"
[ "$1" = "v$CITED" ] || die "tag $1 names a version other than $CITED, which CITATION.cff holds"

printf '%s: %s matches include/version.h and CITATION.cff\n' "$PROGRAM" "$1"
