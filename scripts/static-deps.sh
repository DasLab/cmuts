#!/usr/bin/env bash
#
# static-deps.sh -- build static htslib, HDF5 and their compression libraries
# into one prefix, for `make STATIC_PREFIX=<prefix>`.
#
# usage: scripts/static-deps.sh PREFIX [CC]
#
# The prefix receives include/, lib/ holding each library as a static archive,
# and licenses/ holding each library's license for redistribution. Sources are
# pinned by version and checksum. A library already in lib/ is not rebuilt.
# htslib is built without libcurl, so the binary reads local files only. A musl
# compiler builds everything -static, matching how the binary links.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -euo pipefail

PROGRAM="static-deps.sh"

ZLIB_VERSION=1.3.1
ZLIB_SHA=9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23
ZLIB_URL=https://github.com/madler/zlib/releases/download/v$ZLIB_VERSION/zlib-$ZLIB_VERSION.tar.gz

BZIP2_VERSION=1.0.8
BZIP2_SHA=ab5a03176ee106d3f0fa90e381da478ddae405918153cca248e682cd0c4a2269
BZIP2_URL=https://sourceware.org/pub/bzip2/bzip2-$BZIP2_VERSION.tar.gz

XZ_VERSION=5.6.4
XZ_SHA=269e3f2e512cbd3314849982014dc199a7b2148cf5c91cedc6db629acdf5e09b
XZ_URL=https://github.com/tukaani-project/xz/releases/download/v$XZ_VERSION/xz-$XZ_VERSION.tar.gz

LIBDEFLATE_VERSION=1.24
LIBDEFLATE_SHA=ad8d3723d0065c4723ab738be9723f2ff1cb0f1571e8bfcf0301ff9661f475e8
LIBDEFLATE_URL=https://github.com/ebiggers/libdeflate/archive/refs/tags/v$LIBDEFLATE_VERSION.tar.gz

HTSLIB_VERSION=1.22.1
HTSLIB_SHA=3dfa6eeb71db719907fe3ef7c72cb2ec9965b20b58036547c858c89b58c342f7
HTSLIB_URL=https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2

HDF5_VERSION=1.14.6
HDF5_SHA=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b
HDF5_URL=https://github.com/HDFGroup/hdf5/releases/download/hdf5_$HDF5_VERSION/hdf5-$HDF5_VERSION.tar.gz

MUSL_COPYRIGHT_URL=https://git.musl-libc.org/cgit/musl/plain/COPYRIGHT

die() {
    printf '%s: %s\n' "$PROGRAM" "$1" >&2
    exit 1
}

sha256() {
    if command -v sha256sum > /dev/null; then
        sha256sum "$1"
    else
        shasum -a 256 "$1"
    fi | cut -d' ' -f1
}

# Downloads a pinned tarball into src/, keeping one already there. A checksum
# mismatch removes the file and stops the run.
fetch() {
    local url=$1 sum=$2 out=$SRC/$3

    [ -f "$out" ] || curl -fsSL --retry 3 -o "$out" "$url"

    if [ "$(sha256 "$out")" != "$sum" ]; then
        rm -f "$out"
        die "checksum mismatch for $3"
    fi
}

# Unpacks src/$1 into a fresh work/$2, dropping $3 leading path components.
extract() {
    rm -rf "${WORK:?}/$2"
    mkdir -p "$WORK/$2"
    tar -xf "$SRC/$1" -C "$WORK/$2" --strip-components="$3"
}

# Copies a build tree's license file beside the archives it covers.
license() {
    local dir=$1 name=$2 file

    for file in LICENSE COPYING LICENSE.txt; do
        if [ -f "$dir/$file" ]; then
            cp "$dir/$file" "$LICENSES/$name"
            return
        fi
    done

    die "no license file under $dir"
}

# The binary carries musl's code when built -static with a musl compiler, so
# its copyright notice ships alongside the other licenses. The toolchain or
# the distribution usually installs it; the project's git serves it otherwise.
musl_license() {
    local file

    for file in /usr/share/licenses/musl*/COPYRIGHT /usr/share/licenses/musl*/LICENSE \
                "$(dirname "$(command -v "${CC%% *}")")/../COPYRIGHT"; do
        if [ -f "$file" ]; then
            cp "$file" "$LICENSES/musl"
            return
        fi
    done

    curl -fsSL --retry 3 -o "$LICENSES/musl" "$MUSL_COPYRIGHT_URL" ||
        die "no musl COPYRIGHT found locally or online"
}

build_zlib() {
    [ -f "$PREFIX/lib/libz.a" ] && return
    fetch "$ZLIB_URL" "$ZLIB_SHA" zlib.tar.gz
    extract zlib.tar.gz zlib 1
    (cd "$WORK/zlib" && CC="$CC" CFLAGS="-O2 $STATIC" \
        ./configure --prefix="$PREFIX" --static \
        && make -j"$JOBS" && make install)
    license "$WORK/zlib" zlib
}

build_bzip2() {
    [ -f "$PREFIX/lib/libbz2.a" ] && return
    fetch "$BZIP2_URL" "$BZIP2_SHA" bzip2.tar.gz
    extract bzip2.tar.gz bzip2 1
    (cd "$WORK/bzip2" && make -j"$JOBS" libbz2.a CC="$CC" \
        CFLAGS="-O2 $STATIC -D_FILE_OFFSET_BITS=64" \
        && cp libbz2.a "$PREFIX/lib/" && cp bzlib.h "$PREFIX/include/")
    license "$WORK/bzip2" bzip2
}

build_xz() {
    [ -f "$PREFIX/lib/liblzma.a" ] && return
    fetch "$XZ_URL" "$XZ_SHA" xz.tar.gz
    extract xz.tar.gz xz 1
    (cd "$WORK/xz" && CC="$CC" CFLAGS="-O2 $STATIC" LDFLAGS="$STATIC" \
        ./configure --prefix="$PREFIX" \
        --disable-shared --enable-static --disable-xz --disable-xzdec \
        --disable-lzmadec --disable-lzmainfo --disable-scripts --disable-doc \
        && make -j"$JOBS" && make install)
    license "$WORK/xz" xz
}

build_libdeflate() {
    [ -f "$PREFIX/lib/libdeflate.a" ] && return
    fetch "$LIBDEFLATE_URL" "$LIBDEFLATE_SHA" libdeflate.tar.gz
    extract libdeflate.tar.gz libdeflate 1
    (cd "$WORK/libdeflate" && cmake -B build -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$STATIC" \
        -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY \
        -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF \
        && cmake --build build -j"$JOBS" && cmake --install build)
    license "$WORK/libdeflate" libdeflate
}

build_htslib() {
    [ -f "$PREFIX/lib/libhts.a" ] && return
    fetch "$HTSLIB_URL" "$HTSLIB_SHA" htslib.tar.bz2
    extract htslib.tar.bz2 htslib 1
    (cd "$WORK/htslib" && CC="$CC" ./configure --prefix="$PREFIX" \
        --disable-libcurl --disable-plugins --with-libdeflate \
        CFLAGS="-O2 $STATIC" \
        CPPFLAGS="-I$PREFIX/include" LDFLAGS="$STATIC -L$PREFIX/lib" \
        && make -j"$JOBS" libhts.a)
    cp "$WORK/htslib/libhts.a" "$PREFIX/lib/"
    cp -r "$WORK/htslib/htslib" "$PREFIX/include/"
    license "$WORK/htslib" htslib
}

build_hdf5() {
    [ -f "$PREFIX/lib/libhdf5.a" ] && return
    fetch "$HDF5_URL" "$HDF5_SHA" hdf5.tar.gz
    extract hdf5.tar.gz hdf5 2
    (cd "$WORK/hdf5" && CC="$CC" CFLAGS="-O2 $STATIC" LDFLAGS="$STATIC" \
        ./configure --prefix="$PREFIX" \
        --enable-static --disable-shared --disable-hl --disable-tests \
        --disable-tools --enable-build-mode=production \
        --with-zlib="$PREFIX" \
        && make -j"$JOBS" && make install)
    license "$WORK/hdf5" hdf5
}

[ $# -ge 1 ] && [ $# -le 2 ] || die "usage: $PROGRAM PREFIX [CC]"

mkdir -p "$1"
PREFIX=$(cd "$1" && pwd)
CC=${2:-cc}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}
SRC=$PREFIX/src
WORK=$PREFIX/work
LICENSES=$PREFIX/licenses

command -v "${CC%% *}" > /dev/null || die "no compiler named $CC"
command -v cmake > /dev/null || die "cmake is required, for libdeflate"

# A musl toolchain's test programs only run on the host when linked -static,
# and the final binary links -static under musl anyway.
STATIC=
case "$($CC -dumpmachine)" in
    *musl*) STATIC=-static ;;
esac

mkdir -p "$PREFIX/include" "$PREFIX/lib" "$SRC" "$WORK" "$LICENSES"

build_zlib
build_bzip2
build_xz
build_libdeflate
build_htslib
build_hdf5

[ -z "$STATIC" ] || musl_license

rm -rf "$WORK"
printf '%s: archives ready under %s\n' "$PROGRAM" "$PREFIX/lib"
