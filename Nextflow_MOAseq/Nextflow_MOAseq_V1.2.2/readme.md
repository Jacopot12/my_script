Nextflow pipeline parameters:

**samplesheet**: path to the file.csv with the samples (string)
**genome**: path to the genome file in fasta format (string)
**genome_index**: path to the STAR genome index folder (string), default null
**annotation**: path to the annotation file in gtf format (string), default "NO_FILE"
**outdir**: path to the output directory (string)
**rep_merged**: merge the sample in the samplesheet (boolean), default false
**sample_name**: name of the sample if the samples in the samplesheet are merged (string), default null

_rep_merged_ and _sample_name_ parameters work only in Nextflow_MOAseq_V1.1 and later versions.

All other parameter for the softwares used need to be modify from the modules scripts.
