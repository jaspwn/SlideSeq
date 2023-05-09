
# Troubleshooting

 1. [File not found on NEMO](#file-not-found-on-nemo)
 2. [Uncommitted changes error](#uncommitted-changes-error)
 3. [Wrong puck](#wrong-puck)
 4. [Exotic GTF file](#exotic-gtf-file)
 5. [No usable reads](#no-usable-reads)


## File not found on NEMO

In order to simplify usage of absolute paths, it was decided to mount the whole filesytem when using `singularity` ([here](https://github.com/bahnk/SlideSeqFFPE/blob/spatial/conf/process.config#L5)).
The pipeline was developed on CAMP, so if you're on NEMO `singularity` can't find files as the mounting point is wrong.

To fix this, create a configuration file named `singularity.config` with the following content:

```
singularity {
	enabled = true
	runOptions = "-B /nemo --nv"
}
```

Then, pass this configuration file to nextflow when you run the pipeline:

```
module load Nextflow/22.04.0 Singularity/3.6.4
nextflow pull bahnk/SlideSeq -r main
nextflow -c singularity.config run bahnk/SlideSeq -r main -params-file params.yml --design design.csv
```

## Uncommitted changes error

Sometimes when running the pipeline, user can encounter this error:

```
Project `bahnk/SlideSeq` contains uncommitted changes -- Cannot switch to revision: main
```

This error usually happens because `CAMP/NEMO` changes file permissions after `nextflow pull` (I don't know why).
Nextflow detects these changes and prevents the user from running the pipeline.
To fix this, user needs to go into the directory containing nextflow code and tells `git` to ignore these changes:

```
cd $HOME/.nextflow/assets/bahnk/SlideSeq
git config --local core.filemode false
cd -
```

## Wrong puck

If the pipeline matches too few bead barcodes sometimes it is because user specified the wrong puck.
In order to check that the right puck was specified, user can have a look at the [QC histogram  (page 15 of the report)](output.md#fifth-step-barcode-matching-pages-5-13-14-and-15).

## Exotic GTF file

Most of the problems I had with this pipeline was due to an exotic GTF file.
As illustrated [here](https://github.com/bahnk/slideseq-tools/blob/0294b9bc255b44a94d9abe35063e7f0c75119567/src/count.cpp#L123) and [here](https://github.com/bahnk/slideseq-tools/blob/main/lib/gene.cpp#L26) the pipeline expects the mapping between `gene_name` and `gene_id` to be a bijection.
However, some GTF files don't respect that.
So, usually I write a script that reformats the GTF file so it respects this assumption.

## No usable reads

The pipeline applies a succession of filtering steps in order to select usable reads.
If a sample doesn't contain usable reads the pipeline will fail because the file to process will be empty due to filtering.
So, when the pipeline fails at a certain step, it is good practice to check if the FastQ or BAM file to process is not empty.

