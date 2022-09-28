import java.nio.file.Paths

process umis_per_barcode {

	label "slideseq_tools"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/06_umis_per_barcode/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		
		
		name = metadata["name"]
		suffix = "umis"
		basename = "${name}.${suffix}"
		threshold = params.umis_threshold

		"""
		umis_per_barcode --threshold $threshold $bam "${basename}.bam"
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process bam_metrics {

	label "python"
	label "bam_metrics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/06_umis_per_barcode/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		sed -i 's/^/UMIs per barcode,${name},/g' "${basename}.csv"
		"""
}

process bam_filter {

	label "samtools"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/06_umis_per_barcode/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), val(expr)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		samtools view --expr '${expr}' --output "${basename}.bam" $bam
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process plot_1_arg_1_val {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/06_umis_per_barcode/${filename}" }

	input:
		tuple val(metadata), path(arg), val(value), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg $basename $value
		"""
}

