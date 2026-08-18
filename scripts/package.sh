#!/usr/bin/env bash
#
# package.sh -- assemble a release tarball holding the cmuts binary, the
# cmuts-align script, and every license the binary carries.
#
# usage: scripts/package.sh PREFIX [BUILD]
#
# PREFIX is the directory scripts/static-deps.sh built, read here for the
# licenses of the libraries linked into the binary. BUILD is the directory
# holding the binary, build/release when not given. The tarball is named
# cmuts-<version>-<os>-<arch>.tar.gz, with the version read from the header
# every program reports.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -euo pipefail

PROGRAM="package.sh"

die() {
    printf '%s: %s\n' "$PROGRAM" "$1" >&2
    exit 1
}

[ $# -ge 1 ] && [ $# -le 2 ] || die "usage: $PROGRAM PREFIX [BUILD]"

BUILD=${2:-build/release}

[ -d "$1/licenses" ] || die "$1 holds no licenses; run scripts/static-deps.sh first"
[ -x "$BUILD/cmuts" ] || die "no binary under $BUILD; run make STATIC_PREFIX=$1 first"

VERSION=$(sed -n 's/.*CMUTS_VERSION "\(.*\)".*/\1/p' include/version.h)
OS=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)
NAME=cmuts-$VERSION-$OS-$ARCH

rm -rf "$NAME"
mkdir "$NAME"
cp "$BUILD/cmuts" "$BUILD/cmuts-align" LICENSE "$NAME/"

for file in "$1"/licenses/*; do
    printf '======== %s ========\n\n' "$(basename "$file")"
    cat "$file"
    printf '\n'
done > "$NAME/THIRD_PARTY_LICENSES"

tar -czf "$NAME.tar.gz" "$NAME"
rm -rf "$NAME"
printf '%s: wrote %s.tar.gz\n' "$PROGRAM" "$NAME"
