import java.nio.file.Paths

process shuffling {

	tag { "${name}" }

	label "python"

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf(".ordered.txt") != -1 )
			{
				"${name}/08_puck_barcodes/${filename}"
			}
			else if ( filename.indexOf(".shuffled.txt") != -1 )
			{
				"${name}/08_puck_barcodes/${filename}"
			}
		}

	input:
		tuple val(metadata), path(seq_barcodes), val(puck_barcodes), path(script)

	output:
		tuple val(metadata), path(seq_barcodes), val("ordered"), path("${name}.ordered.txt"), emit: ordered
		tuple val(metadata), path(seq_barcodes), val("shuffled"), path("${name}.shuffled.txt"), emit: shuffled

	script:		

		name = metadata["name"]

		"""
		python3 $script "${puck_barcodes}" "${name}"
		"""
}

