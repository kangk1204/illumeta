        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        mkdir -p idat
        tail -n +2 raw_idat_filelist.tsv | while IFS=$'	' read -r geo basename rest; do
          grn_url=$(printf '%s
' "$rest" | awk -F'	' '{print $(NF-1)}')
          red_url=$(printf '%s
' "$rest" | awk -F'	' '{print $NF}')
          for url in "$grn_url" "$red_url"; do
            [ -n "$url" ] || continue
            out="idat/$(basename "$url")"
            [ -s "$out" ] || curl -L --fail --retry 3 --retry-delay 5 -o "$out" "$url"
          done
        done
