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

# Your code start from here:
trim_galore \
--phred33 \
--clip_R1 20 --clip_R2 20 \
--fastqc \
--fastqc_args "--noextract --outdir cache/ " \
-o cache/ \
--paired input/B73_chr1_subset_reads_1.fastq input/B73_chr1_subset_reads_2.fastq