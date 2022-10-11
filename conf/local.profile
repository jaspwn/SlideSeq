
/*
 * configuration for local execution
 */

profiles {

	local {

		includeConfig "singularity.config"
		
		process.executor = "local"
		
		// processes default configuration
		includeConfig "process/resources.config"
		includeConfig "process/error.config"
		includeConfig "process/publish.config"
		includeConfig "process/containers.config"
		includeConfig "process/gpu.config"
		includeConfig "process/specific.config"
	}
}

