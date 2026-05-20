#!/usr/bin/env bash
# Headless smoke for MetaGap. Usage: bash scripts/smoke.sh [URL]
set -euo pipefail

BASE=${1:-http://127.0.0.1:8000}
COOKIE=$(mktemp)
trap 'rm -f $COOKIE' EXIT

curl -fsSL "$BASE/" -o /dev/null
echo "OK   homepage 200"
curl -fsSL "$BASE/en/login/" -o /dev/null
echo "OK   login page 200"
curl -c "$COOKIE" -b "$COOKIE" -fsSL -X POST "$BASE/en/login/" \
  -d "username=demo_lab&password=metagap-demo" -o /dev/null
echo "OK   login POST"
curl -b "$COOKIE" -fsSL "$BASE/" -o /dev/null
echo "OK   home authed"

echo "SMOKE PASS"
