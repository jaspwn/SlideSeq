
# Running the pipeline

If you're from the Crick, just `ssh` to CAMP, then configure your `params.yml` file and `design.csv` file as detailed [here](config.md), and finally run:

```
module load Nextflow/22.04.0 Singularity/3.6.4
nextflow pull bahnk/SlideSeq -r main
nextflow run bahnk/SlideSeq -r main -params-file params.yml --design design.csv
```

## If you're on NEMO

In order to simplify usage of absolute paths, it was decided to mount the whole filesytem when using `singularity` ([here](https://github.com/bahnk/SlideSeqFFPE/blob/spatial/conf/process.config#L5)).
The pipeline was developed on CAMP, so if you're on NEMO you should use the following procedure.

First, create a configuration file named `singularity.config` with the following content:

```
singularity {
	enabled = true
	runOptions = "-B /nemo --nv"
}
```

Second, pass this configuration file to nextflow when you run the pipeline:

```
module load Nextflow/22.04.0 Singularity/3.6.4
nextflow pull bahnk/SlideSeq -r main
nextflow -c singularity.config run bahnk/SlideSeq -r main -params-file params.yml --design design.csv
```

