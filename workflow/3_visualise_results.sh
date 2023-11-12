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

BG2BW=$1

# Required awk functions
awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $8/$9*100, $5
}
'

awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $4 > "cache/BSMAP_out_"$5".bedGraph"
}
'

#pipe bedGraph to split by context (use dash to read from sdtin)
# per context
awk -F$"\\t" "$awk_make_bedGraph" "cache/BSMAP_out.txt" | \
gawk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

# Make bigWigs files for IGV per context
$BG2BW/bedGraphToBigWig cache/BSMAP_out_CG.bedGraph input/maize_chr1_reference.chrom.sizes output/BSMAP_out_CG.bigWig

$BG2BW/bedGraphToBigWig cache/BSMAP_out_CHG.bedGraph input/maize_chr1_reference.chrom.sizes output/BSMAP_out_CHG.bigWig

$BG2BW/bedGraphToBigWig cache/BSMAP_out_CHH.bedGraph input/maize_chr1_reference.chrom.sizes output/BSMAP_out_CHH.bigWig

# remove intermediate files      
echo "Removing:"
rm -rv cache/BSMAP_out*.bedGraph
