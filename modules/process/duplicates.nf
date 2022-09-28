import java.nio.file.Paths

process select {

	label "slideseq_tools"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/12_duplicates/${filename}" }

	input:
		tuple val(metadata), path(bam), val(suffix)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai"), emit: bam
		tuple val(metadata), val("unique"), path("${basename}.unique.csv"), emit: unique_reads
		tuple val(metadata), val("resolved"), path("${basename}.resolved.csv"), emit: resolved_reads
		tuple val(metadata), val("unresolved"), path("${basename}.unresolved.csv"), emit: unresolved_reads

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"


		"""
		/usr/bin/select "${basename}" $bam "${basename}.bam"

		sed -i '1{s/^/Process,Sample,/}' "${basename}.unique.csv"
		sed -i '2,\${s/^/Duplicates selection,${name},/}' "${basename}.unique.csv"

		sed -i '1{s/^/Process,Sample,/}' "${basename}.resolved.csv"
		sed -i '2,\${s/^/Duplicates selection,${name},/}' "${basename}.resolved.csv"

		sed -i '1{s/^/Process,Sample,/}' "${basename}.unresolved.csv"
		sed -i '2,\${s/^/Duplicates selection,${name},/}' "${basename}.unresolved.csv"

		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process duplicates {

	label "python"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/12_duplicates/${filename}" }

	input:
		tuple val(metadata), path(csv), path(fastq1), path(fastq2), path(script)

	output:
		tuple val(metadata), file("${basename}.csv")

	script:		

		name = metadata["name"]
		status = metadata["status"]
		basename = "${name}.${status}"
		
		"""
		python3 $script $csv $fastq1 $fastq2 "${basename}.csv"
		sed -i '1{s/^/Process,Sample,/}' "${basename}.csv"
		sed -i '2,\${s/^/Duplicates selection,${name},/}' "${basename}.csv"
		"""
}

process bam_metrics {

	label "python"
	label "bam_metrics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/12_duplicates/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		sed -i 's/^/Duplicates selection,${name},/g' "${basename}.csv"
		"""
}

process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/12_duplicates/${filename}" }

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
		saveAs: { filename -> "${name}/12_duplicates/${filename}" }

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

