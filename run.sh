#!/bin/bash -eu

BUILD_DIR="build"

if [[ -d "$BUILD_DIR" ]]; then
    rm -rf "$BUILD_DIR"
fi

mkdir "$BUILD_DIR"
cd "$BUILD_DIR"
cmake ..
make -j
cd test
cp ../../test/data.txt .
./test