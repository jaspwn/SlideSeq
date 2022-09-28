import java.nio.file.Paths

process fastqc {

	label "sequencing"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/00_fastqc/${filename}" }

	input:
		tuple val(metadata), path(fastq1), path(fastq2)

	output:
		tuple val(metadata), file("*_fastqc.html"), emit: html
		tuple val(metadata), file("*_fastqc.zip"), emit: zip

	script:		

		name = metadata["name"]
		
		"""
		fastqc $fastq1 $fastq2
		"""
}

process mark_duplicates {

	label "sequencing"
	time "10:00:00"

	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/04_mark_dulicates/${filename}" }

	input:
		tuple val(metadata), path(bam)
	
	output:
		tuple val(metadata), path("${name}.dup.bam"), emit: bam
		tuple val(metadata), path("${name}.dup.txt"), emit: metrics
	
	script:
		
		name = metadata["name"]

		"""
		picard-tools \
			MarkDuplicates \
				I=$bam \
				O=${name}.dup.bam \
				M=${name}.dup.txt
		"""
}

