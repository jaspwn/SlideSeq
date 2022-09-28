import java.nio.file.Paths

process dge {

	label "slideseq_tools"

	tag { "${name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(bam)

	output:
		tuple val(metadata), file("${directory}")

	script:		
		
		name = metadata["name"]
		gtf = metadata["gtf"]
		directory = "${name}"

		"""
		count --directory . $gtf $bam

		mkdir $directory
		cat matrix.mtx | gzip -c > $directory/matrix.mtx.gz
		cat features.tsv | gzip -c > $directory/features.tsv.gz
		cat barcodes.tsv | gzip -c > $directory/barcodes.tsv.gz
		"""
}

