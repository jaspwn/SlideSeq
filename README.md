
# Slide-Seq

This pipeline pre-configured for people working at the [Francis Crick Institue](https://www.crick.ac.uk/?gclid=EAIaIQobChMIodDA66K59wIVF-vtCh3_SwEJEAAYAiAAEgKrkvD_BwE).
If you're not from the Crick, it should be easy to configure the pipeline for your system.
Just modify these two files:

 * `nexflow.config`
 * `conf/process.config`

Here is the documentation:

 1. [Pipeline structure](doc/structure.md)
 2. [Pipeline configuration](doc/config.md)
 3. [Running the pipeline](doc/run.md)
 4. [Pipeline output](doc/output.md)
 5. [Troubleshooting](doc/troubleshooting.md)


In a nutshell, first be sure that your singularity config directory is not in your home.
For example:

```bash
$ ls -l ~/.singularity
lrwxrwxrwx 1 username domain users 40 Aug 26  2021 /camp/home/username/.singularity -> /camp/stp/babs/working/username/.singularity
```

Then, you can just run:

```bash
# load nextflow and singularity
module load Nextflow/22.04.0 Singularity/3.6.4

# pull the latest version
nextflow pull bahnk/SlideSeq -r main

# download example files
wget https://raw.githubusercontent.com/bahnk/SlideSeq/main/test/design.csv
wget https://raw.githubusercontent.com/bahnk/SlideSeq/main/params.yml

# run the pipeline and pray
nextflow run bahnk/SlideSeq -r main -params-file params.yml --design design.csv
```

If it fails it's probably because your `MODULEPATH` is missing some locations, and/or because you don't have access to the BABS reference area (`/camp/svc/reference/Genomics/babs`).

If you're on NEMO, you need to create a configuration file named `singularity.config`with the following content:

```
singularity {
	enabled = true
	runOptions = "-B /nemo --nv"
}
```

Then, pass this configuration file to nextflow when you run the pipeline:

```
nextflow -c singularity.config run bahnk/SlideSeq -r main -params-file params.yml --design design.csv
```

Any issues => `nourdinebah@gmail.com`. Cheers,

## Comments

The decision was made not to perform deduplication with UMI-tools for now.
Maybe this feature will be added later.
If necessary, UMI-tools can be run afterwards anyway.

