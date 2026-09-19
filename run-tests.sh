#!/bin/sh

set -e

for bin in "$@"
do
    echo "Running ./$bin"
    "./$bin" || exit 1
done
