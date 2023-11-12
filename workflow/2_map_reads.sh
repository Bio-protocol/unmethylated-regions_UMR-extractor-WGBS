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
SAMTOOLS_PATH=$1

# Your code start from here:
# Mapping with BSMAP

bsmap \
-a cache/B73_chr1_subset_reads_1_val_1.fq \
-b cache/B73_chr1_subset_reads_2_val_2.fq \
-d input/maize_chr1_reference.fa \
-o cache/mapped.sam \
-v 5 \
-r 0 \
-q 20

# Sort BSMAP

samtools view -bS cache/mapped.sam > cache/mapped.bam
samtools sort -n cache/mapped.bam -o cache/mapped_nameSrt.bam
samtools fixmate cache/mapped_nameSrt.bam cache/mapped_nameSrt_fixed.bam
samtools sort cache/mapped_nameSrt_fixed.bam -o cache/mapped_sorted.bam
samtools index cache/mapped_sorted.bam

rm cache/mapped.sam cache/mapped_nameSrt.bam \
cache/mapped_nameSrt_fixed.bam cache/mapped.bam

samtools stats cache/mapped_sorted.bam | \
grep ^SN | \
cut -f 2- > output/mapped_stats.txt

bamtools filter \
-isMapped true \
-isPaired true \
-isProperPair true \
-in cache/mapped_sorted.bam \
-out cache/mapped_sorted_pairs.bam

# Remove duplicates
picard MarkDuplicates \
I=cache/mapped_sorted_pairs.bam \
O=cache/mapped_sorted_MarkDup_pairs.bam \
METRICS_FILE=cache/mapped_MarkDupMetrics.txt \
ASSUME_SORTED=true \
CREATE_INDEX=False \
REMOVE_DUPLICATES=true

# Remove overlaps
bam clipOverlap \
--in cache/mapped_sorted_MarkDup_pairs.bam \
--out cache/mapped_sorted_MarkDup_pairs_clipOverlap.bam \
--stats

# cache BSMAP results
python2 lib/methratio.py \
-o cache/methratio.txt \
-d input/maize_chr1_reference.fa \
-u -z \
-s $SAMTOOLS_PATH \
-r cache/mapped_sorted_MarkDup_pairs_clipOverlap.bam

awk -F$"\\t" -f lib/make_bed.awk \
        "cache/methratio.txt" > "cache/BSMAP_out.txt"
