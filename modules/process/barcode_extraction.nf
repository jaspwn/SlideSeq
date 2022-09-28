import java.nio.file.Paths

process extract_barcode {

	label "slideseq_tools"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/02_barcode_extraction/${filename}" }

	input:
		tuple val(metadata), path(read1), path(read2)

	output:
		tuple val(metadata), path("${name}.fastq.gz"), emit: fastq
		tuple val(metadata), path("${name}.extract_barcode.csv"), emit: metrics
		tuple val(metadata), path("${name}.extract_barcode.up_distances.csv"), emit: distances

	script:		
		
		name = metadata["name"]
		read_structure = metadata["read_structure"]

		"""
		extract_barcode \
			--sample "${name}" \
			--fastq "${name}.fastq.gz" \
			--read-structure "${read_structure}" \
			--max-distance ${params.up_errors_threshold} \
			$read1 $read2

		sed -i 's/^/Barcode extraction,/g' "${name}.extract_barcode.csv"
		sed -i 's/^/Barcode extraction,/g' "${name}.extract_barcode.up_distances.csv"
		"""
}

process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/02_barcode_extraction/${filename}" }

	input:
		tuple val(metadata), path(arg), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg $basename
		"""
}

process plot_1_arg_1_value {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/02_barcode_extraction/${filename}" }

	input:
		tuple val(metadata), path(arg), val(value), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg $value $basename
		"""
}

process pcr_duplicates {

	label "slideseq_tools"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/02_barcode_extraction/${filename}" }

	input:
		tuple val(metadata), path(fastq)

	output:
		tuple val(metadata), path("${name}.remove_dups.fastq.gz"), emit: fastq
		tuple val(metadata), path("${name}.remove_dups.csv"), emit: csv

	script:		
		
		name = metadata["name"]

		"""
		pcr_duplicates --sample "${name}" $fastq
		sed -i 's/^/PCR duplicates,/g' "${name}.remove_dups.csv"
		"""
}

