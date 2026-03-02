#!/bin/bash
set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "$0")/../../" && pwd)/tools"
mkdir -p "$TOOLS_DIR"

if [ -f "$TOOLS_DIR/germline/bin/germline" ] || [ -f "$TOOLS_DIR/germline/germline" ]; then
  echo 'GERMLINE already installed'
  exit 0
fi

cd "$TOOLS_DIR"
# wget https://github.com/gusevlab/germline/archive/refs/heads/master.zip  # didn't work on macOS
git clone https://github.com/gusevlab/germline.git 2>/dev/null || true

cd ${TOOLS_DIR}/germline
make clean 2>/dev/null || true
make all

# TODO maybe check for g++ first
if [ -f "./germline" ]; then
  echo "Compiled: $TOOLS_DIR/germline/germline"
  ./germline --help 2>&1 | head -20 || true
elif [ -f "./bin/germline" ]; then
  echo "Compiled: ${TOOLS_DIR}/germline/bin/germline"
  ./bin/germline --help 2>&1 | head -20 || true
else
  echo "ERROR: compilation failed"
  exit 1
fi
