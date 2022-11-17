
# Pipeline structure

The main processes of the pipeline are the following:

 1. We merge the `FASTQ` files associated with the same samples (`merge_lanes`)
 2. We run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (`fastqc`)
 3. We extract the bead barcodes, UP primers and UMIs from read 1 and add this information in the read 2 `FASTQ` file in order to have a unique file to align (`extract_barcode`)
 4. We remove PCR duplicates that same bead barcode, same UMI and same read 2 sequence (`pcr_duplicates`)
 5. We align the `FASTQ` previously generated with [STAR](https://github.com/alexdobin/STAR) (`star`)
 6. We add the bead barcode, UP primer and UMI information as tags in the `BAM` file generated with STAR (`tag_bam`)
 7. We remove the reads that don't match a puck barcode or don't map to the genome (`bam_filter_up_matched`)
 8. We remove the reads whose barcode have less than a given number of UMIs (`bam_filter_umi_threshold`)
 9. For each sequencing barcode, we find the puck barcodes with the smallest Hamming distance (`hamming`)
 10. We try to match the sequencing barcodes with their counterpart in the puck (`matcher`)
 11. We remove the reads that don't match a puck barcode (`bam_filter_barcode_matched`)
 12. We annotate the reads with a gene using [HTSeq](https://htseq.readthedocs.io/en/master/) (`htseq`)
 13. We remove the reads that can't be annotated with a gene (`bam_filter_gene_tags`)
 14. We remove the position duplicates (same bead barcode, same UMI and same mapping position) (`position_duplicates`)
 15. We perform a deduplication step which selects a unique mapping for each (bead barcode,UMI) pair (`select`)
 16. We create a digital expression matrix (`dge`)

