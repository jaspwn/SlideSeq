import java.nio.file.Paths
import java.text.SimpleDateFormat

/////////////////////////////////
def addValue(map, key, value) {//
/////////////////////////////////

	def new_map = map.clone()

	new_map.put(key, value)

	return new_map
}

/////////////////////////////
def removeKeys(map, keys) {//
/////////////////////////////

	def new_map = [:]

	map.each{
		if ( ! keys.contains(it.key) )
		{
			new_map.put(it.key, it.value)
		}
	}

	return new_map
}

///////////////////////////////
def getMinLength(structure) {//
///////////////////////////////

	return structure.split("[A-Z]").collect{it as int }.sum()
}

/////////////////////////
def getPuckName(item) {//
/////////////////////////

	// Here, it's a sun.nio.fs.UnixPath
	def f = new File( item["puck_path"].toString() )

	return f.getName().toString().replaceAll('\\.csv$', '') // single quotes!
}

/////////////////////////////////////
def getDataPaths(item, project_dir)//
/////////////////////////////////////
{
	def map = [:]

	item.each{

		if ( ["fastq_1", "fastq_2", "puck"].contains(it.key) )
		{
			if ( it.value.startsWith("http") )
			{
				// create tmp directory
				def user = System.getenv("USER") ? System.getenv("USER") : "slidesequser"
				def date = new Date()
				def sdf = new SimpleDateFormat("yyyyMMdd")
				def date_str = sdf.format(date)
				def dir_path = Paths.get(project_dir, "${user}-slideseq-${date_str}")
				def dir = new File(dir_path.toString())
				if ( ! dir.exists() ) { dir.mkdirs() }

				// download file if necessary
				def filename = it.value.split("/").last()
				def f = new File( Paths.get(dir.toString(), filename).toString() )
				if ( ! f.exists() )
				{
					println "Downloading ${f}"
					f.bytes = new byte[0]
					def url_matrix = new URL(it.value)
					f << url_matrix.getBytes()
				}
				else
				{
					//println "${f} already here"
				}

				// add path to metadata
				map.put(it.key, Paths.get(f.getAbsolutePath()))
			}

			else
			{
				def f = new File(it.value)
				def full_path = f.getAbsolutePath()
				map.put(it.key, Paths.get(full_path))
			}
		}

		else
		{
			map.put(it.key, it.value)
		}
	}

	map["puck_path"] = map["puck"]
	map.remove("puck")

	return map
}

