#!/bin/sh

####################################
# normalize your dir to the root level of the repo
function canonicalPath
{
    local path="$1" ; shift
    if [ -d "$path" ]
    then
        echo "$(cd "$path" ; pwd)"
    else
        local b=$(basename "$path")
        local p=$(dirname "$path")
        echo "$(cd "$p" ; pwd)/$b"
    fi
}
####################################
# usage example
mycanonicalpath=$(canonicalPath "$../")


# summarise into 100bp tiles
perl lib/met_context_window.pl cache/BSMAP_out.txt 100

R -f lib/Call-umrs.R \
--args input/maize_chr1_reference_100bp_tiles_sites_counts.txt \
3 \
2 \
0.4 \
0.1 \
cache \
output
