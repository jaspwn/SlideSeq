
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

	def f = new File(puck)

	return f.getName().toString().replaceAll('\\.csv$', '') // single quotes!
}


