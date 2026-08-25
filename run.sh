#!/usr/bin/env sh
set -eu
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "${JOBS:-2}"
exec python3 tools/server.py "$@"
