
# Running the pipeline

If you're from the Crick, just `ssh` to CAMP, then configure your `params.yml` file and `design.csv` file as detailed [here](doc/config.md), and finally run:

```
module load Nextflow/22.04.0 Singularity/3.6.4
nextflow pull bahnk/SlideSeq -r main
nextflow run bahnk/SlideSeq -r main -params-file params.yml --sample_sheet design.csv
```
