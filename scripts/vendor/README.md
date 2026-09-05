# WES annotation recovery utility

`vcf2maf-1.6.22.pl` is vendored without modification from [mskcc/vcf2maf v1.6.22](https://github.com/mskcc/vcf2maf/tree/v1.6.22); see `vcf2maf-LICENSE.txt` (Apache 2.0). Downloaded 2026-09-05.

The annotated TOV3121D VCF was present in the raw archive but its MAF conversion was missing. `scripts/lib/wes_maf_inputs.R` invokes this converter with `--inhibit-vep`, retaining original Ensembl VEP113 annotation and using GRCh38. No variant calling or consequence annotation is repeated. `scripts/lib/twobit_faidx.R` implements the reference-extraction operation used by the converter using the installed hg38 TwoBit genome and rtracklayer; it is not a general samtools replacement. Tabix is not invoked in this mode. The source MAF/VCF paths and SHA256 hashes are recorded in `output/wes_input_manifest.csv`.

Validation against the previously available TOV2835EP MAF reproduced all29,323 records exactly across gene/coordinates/alleles/classification/consequence/protein-position/population-AF/ClinVar fields. The recovered MAF uses the same historical tumour-barcode convention; quantitative VAF is read directly from the archived VCF in script07.
