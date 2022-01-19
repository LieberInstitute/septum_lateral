## Notes Re: checking R1 read lengths:

Usually with FASTQ files from the JHU Single Cell & Transcriptomics Core, it's recommended each `R1` `.fastq` file gets checked for all reads being a fixed length of 28bp.  This seemed to not have been trimmed as with data from the core before, so memory ran out in tabulating these read lengths, in a `qsub` batch `R` job.

However, since `cellranger count` allows you to hard-trim read inputs, just used this utility with `--r1-length=28`, and the resulting Cell Ranger output looked good.

This 'step' can thus be ignored. -MNT 19Jan2022
