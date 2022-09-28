import java.nio.file.Paths

process hamming {

	label "gpu"
	
	time "06:00:00"
	memory "40G"

	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/09_hamming_distance/${filename}" }

	input:
		tuple val(metadata), path(read_barcodes), path(puck_barcodes)

	output:
		tuple val(metadata), file("${csv}")

	script:		
		
		name = metadata["name"]
		bcd = metadata["barcodes"]
		basename = "${name}.${bcd}"
		csv = "${basename}.hamming.csv"

		"""
		hostname
		hamming $read_barcodes $puck_barcodes "${csv}"
		sed -i 's/^/Hamming distance,${name},/g' "${csv}"
		"""
}

process plot_2_args {

	label "plot"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "temp_files" ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename -> "${name}/09_hamming_distance/${filename}" }

	input:
		tuple val(metadata), path(arg1), path(arg2), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg1 $arg2 $basename
		"""
}

