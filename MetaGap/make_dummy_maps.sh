#!/usr/bin/env bash
set -euo pipefail

# root staticfiles path (adjust if needed)
STATIC_ROOT="staticfiles/app"

echo "Creating dummy .map files for CSS/JS in $STATIC_ROOT …"

find "$STATIC_ROOT" \( -name "*.css" -o -name "*.js" \) | while read -r file; do
  dir=$(dirname "$file")
  base=$(basename "$file")

  # initialize mapfile vars
  mapfile1=""
  mapfile2=""
  mapfile3=""

  # always possible fallback: file + ".map"
  mapfile1="${file}.map"

  # CSS cases
  if [[ "$base" == *.min.css ]]; then
    mapfile2="${dir}/${base}.map"
  elif [[ "$base" == *.css ]]; then
    mapfile2="${dir}/${base}.map"
  fi

  # JS cases
  if [[ "$base" == *.min.js ]]; then
    mapfile3="${dir}/${base}.map"
  elif [[ "$base" == *.js ]]; then
    mapfile3="${dir}/${base}.map"
  fi

  for m in "$mapfile1" "$mapfile2" "$mapfile3"; do
    # skip empty values
    if [ -n "$m" ] && ! [ -f "$m" ]; then
      mkdir -p "$(dirname "$m")"
      echo "{}" > "$m"   # or leave empty
      echo "Created dummy map: $m"
    fi
  done
done
