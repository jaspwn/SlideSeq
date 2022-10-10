import java.nio.file.Paths

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
def getPuckName(puck) {//
/////////////////////////

	// I don't know if it's a string or a path
	// Here, it's a sun.nio.fs.UnixPath
	def f = new File( puck.toString() )

	return f.getName().toString().replaceAll('\\.csv$', '') // single quotes!
}

/////////////////////////////////
def convertToAbsolutePath(path)//
/////////////////////////////////
{
	def f = new File(path)
	def full_path = f.getAbsolutePath() 
	return Paths.get(full_path)
}

