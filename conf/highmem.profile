
/*
 * configuration for large files with high memory nodes
 */

profiles {

	highmem {

		includeConfig "singularity.config"
		
		process.executor = "slurm"
		
		// processes default configuration
		includeConfig "process/resources.config"
		includeConfig "process/error.config"
		includeConfig "process/publish.config"
		includeConfig "process/containers.config"
		includeConfig "process/gpu.config"
		includeConfig "process/specific.config"
		
		// high memory nodes for certain process
		includeConfig "process/hmem.config"
	}
}

