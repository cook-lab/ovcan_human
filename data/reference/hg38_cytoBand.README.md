# hg38 cytoband reference

Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
Downloaded and decompressed 2026-09-05. SHA256 of `hg38_cytoBand.tsv`: `ce1b6033a5243e7c5022660b952d2ec33243e307e909afcaeec1894641a5208f`.

Columns: chromosome, 0-based start, half-open end, band, stain. Script08 excludes each chromosome's full `acen` interval from both arms and intersects segment lengths with each non-centromeric arm. Acrocentric p arms13/14/15/21/22 are unassessed.
