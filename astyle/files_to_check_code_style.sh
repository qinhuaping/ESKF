#!/usr/bin/env bash
set -eu

PATTERN="-e ."

if [ $# -gt 0 ]
then
    PATTERN="$1"
fi

exec find src \
    -path src/ESKF -prune -o \
    -path src/Node -prune -o \
    -path include/eskf/ESKF -prune -o \
    -path include/eskf/Node -prune -o \
    -type f \( -name "*.c" -o -name "*.h" -o -name "*.cpp" -o -name "*.hpp" \) | grep $PATTERN
