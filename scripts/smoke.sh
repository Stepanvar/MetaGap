#!/usr/bin/env bash
# Headless smoke for MetaGap. Usage: bash scripts/smoke.sh [URL]
set -euo pipefail

BASE=${1:-http://127.0.0.1:8000}
COOKIE=$(mktemp)
trap 'rm -f $COOKIE' EXIT

curl -fsSL "$BASE/" -o /dev/null
echo "OK   homepage 200"
# GET login page first to obtain CSRF cookie
curl -c "$COOKIE" -b "$COOKIE" -fsSL "$BASE/en/login/" -o /dev/null
echo "OK   login page 200"
CSRF=$(grep csrftoken "$COOKIE" | awk '{print $NF}')
# POST without auto-following redirect so session cookie is saved before next request
curl -c "$COOKIE" -b "$COOKIE" -sS -X POST "$BASE/en/login/" \
  -d "username=demo_lab&password=metagap-demo&csrfmiddlewaretoken=$CSRF" \
  -H "Referer: $BASE/en/login/" \
  -o /dev/null -w "%{http_code}\n" | grep -qE "^30[12]$"
echo "OK   login POST"
curl -b "$COOKIE" -fsSL "$BASE/" -o /dev/null
echo "OK   home authed"

echo "SMOKE PASS"
