import java.nio.file.Paths

process matcher {

	label "python"
	time "02:00:00"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/10_barcode_matching/${filename}" }

	input:
		tuple val(metadata), path(hamming), path(reads), path(puck), path(script)

	output:
		tuple val(metadata), file("${basename}.map.matching.csv"), emit: mapping
		tuple val(metadata), file("${basename}.values.matching.csv"), emit: values
		tuple val(metadata), file("${basename}.metrics.matching.csv"), emit: metrics
		tuple val(metadata), file("${basename}.csv"), emit: coords

	script:		
		
		name = metadata["name"]
		shuf = metadata["barcodes"]
		basename = name + "." + shuf

		"""
		python3 $script $hamming $reads $puck "${basename}" \
			$params.barcode_errors_threshold \
			$params.barcode_max_matches \
			$params.barcode_max_entropy
		"""
}

process add_match {

	label "slideseq_tools"
	time "03:00:00"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/10_barcode_matching/${filename}" }

	input:
		tuple val(metadata), path(mapping), path(bam)

	output:
		tuple \
			val(metadata),
			file("${name}.matched.bam"),
			file("${name}.matched.bam.bai")

	script:		
		
		name = metadata["name"]

		// add_match requires the csv file to have only 2 columns
		// so we remove the first and second ones

		"""
		cut -d "," -f 3,4 $mapping > tmp.csv
		add_match tmp.csv $bam "${name}.matched.bam"
		samtools index "${name}.matched.bam"
		"""
}

process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/10_barcode_matching/${filename}" }

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

process bam_metrics {

	label "python"
	label "bam_metrics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/10_barcode_matching/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		sed -i 's/^/Barcode matching,${name},/g' "${basename}.csv"
		"""
}

process bam_filter {

	label "samtools"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/10_barcode_matching/${filename}" }

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

