import java.nio.file.Paths

process get_barcodes {

	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/07_sequencing_barcodes/${filename}" }

	input:
		tuple val(metadata), path(csv)

	output:
		tuple val(metadata), file("${name}.barcodes.txt")

	script:		
		
		name = metadata["name"]

		"""
		cut -d "," -f 3 $csv | sort | uniq > "${name}.barcodes.txt"
		"""
}

