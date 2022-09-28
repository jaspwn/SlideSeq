#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths

///////////////////////////////////////////////////////////////////////////////
//// METHODS //////////////////////////////////////////////////////////////////

include {
	addValue;
	removeKeys;
	getMinLength;
	getPuckName
	} from "./modules/utils.nf"

///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////

//////////////////
// quality control
include { fastqc } from "./modules/process/quality_control"
include { mark_duplicates } from "./modules/process/quality_control"
//////////////////

//////////
// samples
include { merge_lanes } from "./modules/process/preprocessing"
//////////

/////////////////////
// barcode extraction

include { extract_barcode } from "./modules/process/barcode_extraction"

include { plot_1_arg as plot_up_matching } from "./modules/process/barcode_extraction"
plot_up_matching_script = Channel.fromPath("bin/plot/up_matching.py")

include { plot_1_arg_1_value as plot_barcode_extraction } from "./modules/process/barcode_extraction"
plot_barcode_extraction_script  = Channel.fromPath("bin/plot/barcode_extraction.py")
/////////////////////

///////////////////////////
// alignment and duplicates
include { star } from "./modules/process/alignment"
///////////////////////////

/////////////////
// slide seq tags
include { bam_tag as tag_bam } from "./modules/process/bam_tags"

include { bam_metrics as reads_up_matching } from "./modules/process/bam_tags"
reads_up_matching_script = Channel.fromPath("bin/bam/reads_up_matching.py")

include { plot_1_arg as plot_up_align } from "./modules/process/bam_tags"
plot_up_align_script = Channel.fromPath("bin/plot/up_align.py")

include { bam_filter as bam_filter_up_matched } from "./modules/process/bam_tags"
/////////////////

/////////////////////////////
// umis per barcode threshold

include { umis_per_barcode } from "./modules/process/umis_per_barcode"

include { bam_metrics as reads_umis_per_barcode } from "./modules/process/umis_per_barcode"
reads_umis_per_barcode_script = Channel.fromPath("bin/bam/reads_umis_per_barcode.py")

include { bam_metrics as reads_umi_threshold } from "./modules/process/umis_per_barcode"
reads_umi_threshold_script = Channel.fromPath("bin/bam/reads_umi_threshold.py")

include { bam_filter as bam_filter_umi_threshold } from "./modules/process/umis_per_barcode"

include { plot_1_arg_1_val as plot_umi_threshold } from "./modules/process/umis_per_barcode"
plot_umi_threshold_script = Channel.fromPath("bin/plot/umi_threshold.py")
/////////////////////////////

////////////////
// puck barcodes
include { shuffling } from "./modules/process/puck_barcodes"
shuffling_script = Channel.fromPath("bin/shuffling.py")
////////////////

//////////////////////
// sequencing barcodes
include { get_barcodes } from "./modules/process/sequencing_barcodes"
//////////////////////

///////////////////
// hamming distance

include { hamming } from "./modules/process/hamming"

include { plot_2_args as plot_histo_hamming } from "./modules/process/hamming"
plot_histo_hamming_script = Channel.fromPath("bin/plot/histo_hamming.py")
///////////////////

///////////////////
// barcode matching

include { matcher } from "./modules/process/barcode_matching"
matcher_script = Channel.fromPath("bin/matcher.py")

include { add_match } from "./modules/process/barcode_matching"

include { plot_1_arg as plot_barcode_matching } from "./modules/process/barcode_matching"
plot_barcode_matching_script = Channel.fromPath("bin/plot/barcode_matching.py")

include { plot_1_arg as plot_barcode_align } from "./modules/process/barcode_matching"
plot_barcode_align_script = Channel.fromPath("bin/plot/barcode_align.py")

include { plot_1_arg as plot_histo_errors } from "./modules/process/barcode_matching"
plot_histo_errors_script = Channel.fromPath("bin/plot/histo_errors.py")

include { bam_metrics as reads_barcode_matching } from "./modules/process/barcode_matching"
reads_barcode_matching_script = Channel.fromPath("bin/bam/reads_barcode_matching.py")

include { bam_filter as bam_filter_barcode_matched } from "./modules/process/barcode_matching"
///////////////////

///////////////
// gene tagging

include { htseq } from "./modules/process/gene_tags"

include { bam_metrics as count_gene_tags } from "./modules/process/gene_tags"
count_gene_tags_script = Channel.fromPath("bin/bam/count_gene_tags.py")

include { bam_metrics as count_reads_per_umi } from "./modules/process/gene_tags"
count_reads_per_umi_script = Channel.fromPath("bin/bam/reads_per_umi.py")

include { bam_metrics as count_reads_per_umi_gene } from "./modules/process/gene_tags"
count_reads_per_umi_gene_script = Channel.fromPath("bin/bam/reads_per_umi_gene.py")

include { plot_1_arg as plot_gene_tags } from "./modules/process/gene_tags"
plot_gene_tags_script = Channel.fromPath("bin/plot/gene_tags.py")

include { bam_filter as bam_filter_gene_tags } from "./modules/process/gene_tags"
///////////////

/////////////
// duplicates

include { select } from "./modules/process/duplicates"

include { duplicates } from "./modules/process/duplicates"
duplicates_script = Channel.fromPath("bin/duplicates.py")

include { bam_metrics as count_select } from "./modules/process/duplicates"
count_select_script = Channel.fromPath("bin/bam/count_select.py")

include { plot_1_arg as plot_select } from "./modules/process/duplicates"
plot_select_script = Channel.fromPath("bin/plot/select.py")

include { bam_filter as bam_filter_multimapped_umis } from "./modules/process/duplicates"
/////////////

///////////////////
// barcodes metrics

include { bam_metrics as reads_per_barcode_umi } from "./modules/process/barcodes_metrics"
reads_per_barcode_umi_script = Channel.fromPath("bin/bam/reads_per_barcode_umi.py")

include { plot_1_arg as plot_balance_barcode } from "./modules/process/barcodes_metrics"
plot_balance_barcode_script = Channel.fromPath("bin/plot/balance_barcode.py")

include { plot_1_arg as plot_balance_umi } from "./modules/process/barcodes_metrics"
plot_balance_umi_script  = Channel.fromPath("bin/plot/balance_umi.py")

include { plot_1_arg as plot_reads_fraction } from "./modules/process/barcodes_metrics"
plot_reads_fraction_script = Channel.fromPath("bin/plot/reads_fraction.py")
///////////////////

//////
// dge
include { dge } from "./modules/process/digital_expression_matrix"
//////

/////////////////
// counts metrics

plot_histo_genes_script = Channel.fromPath("bin/plot/histo_genes.py")
include { plot_1_arg as plot_histo_genes } from "./modules/process/counts_metrics"

plot_histo_umis_script = Channel.fromPath("bin/plot/histo_umis.py")
include { plot_1_arg as plot_histo_umis } from "./modules/process/counts_metrics"

plot_spatial_umis_script  = Channel.fromPath("bin/plot/spatial_umi.py")
include { plot_1_arg_csv_output as plot_umis_per_barcode } from "./modules/process/counts_metrics"

plot_umis_per_barcode_script = Channel.fromPath("bin/plot/umis_per_barcode.py")
include { plot_2_args as plot_spatial_umis } from "./modules/process/counts_metrics"
/////////////////

/////////
// export

include { merge_plots } from "./modules/process/export"
include { reformat_coords } from "./modules/process/export"
include { export_metrics } from "./modules/process/export"
/////////

///////////////////////////////////////////////////////////////////////////////
//// DESIGN ///////////////////////////////////////////////////////////////////

Channel
	.fromPath(params.design)
	.splitCsv(header: true)
	.map{ addValue(it, "min_length", getMinLength(it["read_structure"])) }
	.map{ addValue(it, "puck_path", new File(it["puck"]).getAbsolutePath()) }
	.map{ addValue(it, "puck", getPuckName(it["puck_path"])) }
	.map{ addValue(it, "fastq_1", new File(it["fastq_1"]).getAbsolutePath()) }
	.map{ addValue(it, "fastq_2", new File(it["fastq_2"]).getAbsolutePath()) }
	.set{ FASTQ }

///////////////////////////////////////////////////////////////////////////////
//// PUCKS ////////////////////////////////////////////////////////////////////

FASTQ
	.map{ [ it["puck"] , it["puck_path"] ] }
	.unique()
	.set{ PUCKS }

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////

workflow {

	///////////////////////////////////////////////////////////////////////////
	// FASTQC

	FASTQ
		.map{[
			removeKeys(it, ["fastq_1", "fastq_2"]),
			it["fastq_1"],
			it["fastq_2"]
		]}
		.set{ TO_FASTQC }

	fastqc(TO_FASTQC)

	///////////////////////////////////////////////////////////////////////////
	// MERGE

	FASTQ
		.map{[
			it["name"],
			removeKeys(it, ["fastq_1", "fastq_2"]),
			it["fastq_1"],
			it["fastq_2"]
		]}
		.groupTuple()
		.map{[
			it[1][0],
			it[2].sort{ new File(it).getName() },
			it[3].sort{ new File(it).getName() }
		]}
		.map{ [ addValue(it[0], "n_fastq_read1", it[1].size()) , *it[1..2] ] }
		.map{ [ addValue(it[0], "n_fastq_read2", it[2].size()) , *it[1..2] ] }
		.set{ TO_MERGE }
	
	merge_lanes(TO_MERGE)

	///////////////////////////////////////////////////////////////////////////
	// BARCODE EXTRACTION

	extract_barcode(merge_lanes.out)

	plot_up_matching(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from("up_matching") )
			.combine(plot_up_matching_script)
	)

	plot_barcode_extraction(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from("barcode_extraction") )
			.combine(plot_barcode_extraction_script)
			.map{ [ it[0] , it[1] , it[0]["min_length"] , it[2] , it[3] ] }
	)

	///////////////////////////////////////////////////////////////////////////
	// ALIGNMENT

	star(extract_barcode.out.fastq)

	///////////////////////////////////////////////////////////////////////////
	// DUPLICATES

	mark_duplicates(star.out.bam)

	///////////////////////////////////////////////////////////////////////////
	// ADD SLIDE-SEQ, ALIGNMENT AND DUPLICATES TAGS

	// slide-seq tags, alignment tag and duplicate tag
	tag_bam(
		mark_duplicates
			.out
			.bam
			.combine( Channel.from("tagged") )
	)

	// 3 columns: Matched, Mapped, Reads
	reads_up_matching(
		tag_bam
			.out
			.combine( Channel.from("reads_up_matching") )
			.combine(reads_up_matching_script)
	)

	plot_up_align(
		reads_up_matching
			.out
			.combine( Channel.from("up_align") )
			.combine(plot_up_align_script)
	)

	bam_filter_up_matched(
		tag_bam
			.out
			.combine( Channel.from("up_matched") )
			.combine( Channel.from("[us]==\"MATCHED\" && [as]==\"MAPPED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// UMIS PER BARCODE THRESHOLD

	umis_per_barcode(bam_filter_up_matched.out)

	reads_umis_per_barcode(
		umis_per_barcode
			.out
			.combine( Channel.from("reads_umis_per_barcode") )
			.combine(reads_umis_per_barcode_script)
	)

	reads_umi_threshold(
		umis_per_barcode
			.out
			.combine( Channel.from("reads_umi_threshold") )
			.combine(reads_umi_threshold_script)
	)

	plot_umi_threshold(
		reads_umi_threshold
			.out
			.combine( Channel.from(params.umis_threshold) )
			.combine( Channel.from("umi_threshold") )
			.combine(plot_umi_threshold_script)
	)

	bam_filter_umi_threshold(
		umis_per_barcode
			.out
			.combine( Channel.from("umi_threshold") )
			.combine( Channel.from("[bt]==\"PASS\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// SEQUENCING BARCODES

	get_barcodes(reads_umis_per_barcode.out)

	///////////////////////////////////////////////////////////////////////////
	// PUCK BARCODES

	get_barcodes
		.out
		.combine( PUCKS )
		.filter{ it[0]["puck"] == it[2] }
		.map{ [ *it[0..1] , it[3] ] }
		.set{ TO_SHUFFLING }

	shuffling( TO_SHUFFLING.combine(shuffling_script) )

	///////////////////////////////////////////////////////////////////////////
	// HAMMING DISTANCE

	shuffling
		.out
		.ordered
		.concat( shuffling.out.shuffled )
		.map{ [ addValue(it[0], "barcodes", it[2]) , it[1] , it[3] ] }
		.set{ TO_HAMMING }

	hamming(TO_HAMMING)

	plot_histo_hamming(
		hamming
			.out
			.map{ [ it[0]["name"] , *it ] }
			.groupTuple()
			.map{
				it[1][0]["barcodes"] == "ordered" ? [it[1][0], *it[2]] : [it[1][1], *it[2].reverse()]
			}
			.combine( Channel.from("histo_hamming") )
			.combine(plot_histo_hamming_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// BARCODE MATCHING

	hamming
		.out
		.combine(reads_umis_per_barcode.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ it[0] , it[1] , it[3] ] }
		.combine( PUCKS )
		.filter{ it[0]["puck"] == it[3] }
		.map{ [ *it[0..2] , it[4] ] }
		.set{ TO_MATCHING }
	
	matcher( TO_MATCHING.combine(matcher_script) )

	matcher
		.out
		.mapping
		.filter{ it[0]["barcodes"] == "ordered" }
		.combine(bam_filter_umi_threshold.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ *it[0..1] , it[3] ] }
		.set{ TO_ADD_MATCH }

	plot_barcode_matching(
		matcher
			.out
			.metrics
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("barcode_matching") )
			.combine(plot_barcode_matching_script)
	)

	plot_histo_errors(
		matcher
			.out
			.mapping
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("histo_errors") )
			.combine(plot_histo_errors_script)
	)

	add_match(TO_ADD_MATCH)

	reads_barcode_matching(
		add_match
			.out
			.combine( Channel.from("reads_barcode_matching") )
			.combine(reads_barcode_matching_script)
	)

	plot_barcode_align(
		reads_barcode_matching
			.out
			.combine( Channel.from("barcode_align") )
			.combine(plot_barcode_align_script)
	)

	bam_filter_barcode_matched(
		add_match
			.out
			.combine( Channel.from("barcode_matched") )
			.combine( Channel.from("[bs]==\"MATCHED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// GENE TAGS

	htseq(bam_filter_barcode_matched.out)

	count_gene_tags(
		htseq
			.out
			.bam
			.combine( Channel.from("gene_tags") )
			.combine(count_gene_tags_script)
	)

	plot_gene_tags(
		count_gene_tags
			.out
			.combine( Channel.from("gene_tags") )
			.combine(plot_gene_tags_script)
	)

	bam_filter_gene_tags(
		htseq
			.out
			.bam
			.combine( Channel.from("gene_tags") )
			.combine( Channel.from("[XF]!~\"^__.+\"") )
	)

	count_reads_per_umi(
		bam_filter_gene_tags
			.out
			.combine( Channel.from("reads_per_umi") )
			.combine(count_reads_per_umi_script)
	)

	count_reads_per_umi_gene(
		bam_filter_gene_tags
			.out
			.combine( Channel.from("reads_per_umi_gene") )
			.combine(count_reads_per_umi_gene_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// DUPLICATES

	select(
		bam_filter_gene_tags
			.out
			.map{ it[0..1] }
			.combine( Channel.from("select") )
	)

	select
		.out
		.unique_reads
		.concat(select.out.resolved_reads)
		.concat(select.out.unresolved_reads)
		.combine(merge_lanes.out)
		.filter{ it[0]["name"] == it[3]["name"] }
		.map{ [ addValue(it[0], "status", it[1]) , it[2] , *it[4..5] ] }
		.set{ TO_DUPLICATES }

	duplicates(TO_DUPLICATES.combine(duplicates_script))

	count_select(
		select
			.out
			.bam
			.combine( Channel.from("count_select") )
			.combine(count_select_script)
	)

	plot_select(
		count_select
			.out
			.combine( Channel.from("selected") )
			.combine(plot_select_script)
	)

	bam_filter_multimapped_umis(
		select
			.out
			.bam
			.combine( Channel.from("multimap_umis") )
			.combine( Channel.from("[cs]==\"UNIQUE\" || [cs]==\"INCLUDED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// BARCODES METRICS

	reads_per_barcode_umi(
		bam_filter_multimapped_umis
			.out
			.combine( Channel.from("reads_per_barcode_umi") )
			.combine(reads_per_barcode_umi_script)
	)

	plot_balance_barcode(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("balance_barcode") )
			.combine(plot_balance_barcode_script)
	)

	plot_balance_umi(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("balance_umi") )
			.combine(plot_balance_umi_script)
	)

	plot_reads_fraction(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("reads_fraction") )
			.combine(plot_reads_fraction_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// DIGITAL EXPRESSION MATRIX

	dge( bam_filter_multimapped_umis.out.map{it[0..1]} )

	///////////////////////////////////////////////////////////////////////////
	// COUNTS METRICS

	plot_histo_umis(
		dge
			.out
			.combine( Channel.from("histo_umis") )
			.combine(plot_histo_umis_script)
	)

	plot_histo_genes(
		dge
			.out
			.combine( Channel.from("histo_genes") )
			.combine(plot_histo_genes_script)
	)

	plot_umis_per_barcode(
		dge
			.out
			.combine( Channel.from("umis_per_barcode") )
			.combine(plot_umis_per_barcode_script)
	)

	plot_spatial_umis(
		dge
			.out
			.combine(matcher.out.coords)
			.filter{ it[2]["barcodes"] == "ordered" }
			.filter{ it[0]["name"] == it[2]["name"] }
			.map{ [ * it[0..1] , it[3] ] }
			.combine( Channel.from("spatial_umis") )
			.combine(plot_spatial_umis_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// EXPORT

	// metrics
	extract_barcode.out.metrics
		.concat(	
			extract_barcode.out.distances,
			reads_up_matching.out,
			//reads_umis_per_barcode.out,
			reads_umi_threshold.out,
			matcher.out.metrics,
			reads_barcode_matching.out,
			count_gene_tags.out,
			//count_reads_per_umi.out,
			//count_reads_per_umi_gene.out,
			count_select.out,
			//reads_per_barcode_umi.out,
			plot_umis_per_barcode.out.csv
		)
		.map{ it[1] }
		.collect()
		.set{ ALL_METRICS }
	
	export_metrics(ALL_METRICS)
	
	// plot
	plot_barcode_extraction.out.pdf
		.concat(
			plot_up_matching.out.pdf,
			plot_up_align.out.pdf,
			plot_umi_threshold.out.pdf,
			plot_barcode_align.out.pdf,
			plot_gene_tags.out.pdf,
			plot_select.out.pdf,
			plot_balance_barcode.out.pdf,
			plot_balance_umi.out.pdf,
			plot_reads_fraction.out.pdf,
			plot_histo_umis.out.pdf,
			plot_histo_genes.out.pdf,
			plot_barcode_matching.out.pdf,
			plot_histo_errors.out.pdf,
			plot_histo_hamming.out.pdf,
			plot_umis_per_barcode.out.pdf,
			plot_spatial_umis.out.pdf
		)
		.map{ [ it[0]["name"] , it[1] ] }
		.groupTuple()
		.set{ TO_MERGE_PLOTS }
	
	merge_plots(TO_MERGE_PLOTS)

	// spatial info
	reformat_coords( matcher.out.coords.filter{ it[0]["barcodes"] == "ordered"} )
}

