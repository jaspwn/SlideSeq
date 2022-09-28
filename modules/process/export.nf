import java.nio.file.Paths

process merge_plots {

	label "python"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(name), path(pdfs)

	output:
		file "${name}.pdf"

	script:		
		"""
		pdfunite $pdfs "${name}.pdf"
		"""
}

process reformat_coords {

	tag { "${name}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(csv)

	output:
		file "${name}.csv"

	script:		
		
		name = metadata["name"]

		// The header is this one:
		// Process,Sample,SeqBarcode,PuckBarcode,x,y,Distance,Number,Reads
		// So, we take the 3rd column (SeqBarcode)

		"""
		cut -d "," -f 3,5,6 $csv > "${name}.csv"
		sed -i 's/SeqBarcode/Barcode/g' "${name}.csv"
		"""
}

process export_metrics {

	label "python"
	
	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true"

	input:
		path csvs

	output:
		file "qc_metrics.csv"

	script:		
		"""
		cat $csvs > "qc_metrics.csv"
		"""
}

