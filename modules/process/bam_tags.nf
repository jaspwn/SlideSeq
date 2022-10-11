import java.nio.file.Paths

process bam_tag {

	label "slideseq_tools"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/05_add_bam_tags/${filename}" }

	input:
		tuple val(metadata), path(bam), val(suffix)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"


		"""
		tag_bam $bam "${basename}.bam"
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
		saveAs: { filename -> "${name}/05_add_bam_tags/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		sed -i 's/^/BAM tagging,${name},/g' "${basename}.csv"
		"""
}

process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/05_add_bam_tags/${filename}" }

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

process bam_filter {

	label "samtools"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/05_add_bam_tags/${filename}" }

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

