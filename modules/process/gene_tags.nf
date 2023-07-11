import java.nio.file.Paths

process htseq {

	label "sequencing"

	tag { "${name}" }
	
	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/11_gene_tags/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai)

	output:
		tuple val(metadata), path("${out_bam}"), path("${out_bam}.bai"), emit: bam
		tuple val(metadata), path("${out_txt}"), emit: txt

	script:		
		
		name = metadata["name"]
		out_sam = "${name}.htseq.sam"
		out_bam = "${name}.htseq.bam"
		out_txt = "${name}.htseq.txt"
		gtf = metadata["gtf"]

		"""
		htseq-count -i gid --format bam --samout sample.sam $bam $gtf > $out_txt
		samtools view -H $bam > $out_sam
		cat sample.sam >> $out_sam
		samtools view -S -b $out_sam > $out_bam
		samtools index $out_bam
		"""
}

process bam_metrics {

	label "python"
	label "bam_metrics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/11_gene_tags/${filename}" }

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		sed -i 's/^/Gene tagging,${name},/g' "${basename}.csv"
		"""
}

process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/11_gene_tags/${filename}" }

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
		saveAs: { filename -> "${name}/11_gene_tags/${filename}" }

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

