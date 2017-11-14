# GATK GermlineCNVCaller
Using GATK-4.beta.5

## Get suitable samples

- Ran with seqcapv3
- on-target coverage >20x for at least 75%

Perform gentli query using [find_good_samples.py](https://github.com/wdecoster/GermlineCNVCaller/blob/master/find_good_samples.py)

## Arguments to GermlineCNVCaller
* --jobType  
either LEARN_AND_CALL or CALL_ONLY
* --input  
Combined raw or GC corrected (but not proportional) read counts table. Can be for a cohort or for a single sample. The format is described in ReadCountCollectionUtils  
* --sexGenotypeTable  
Sample sex genotypes table. The format is described in SexGenotypeTableReader, or a table listing the known sex genotypes of each sample. Only two columns, labeled SAMPLE_NAME and SEX_GENOTYPE, are required. The SEX_GENOTYPE column values must match those of the contigAnnotationsTable.
`cat test_targetcoverage | grep contig | head -n 1 | cut -f4- | tr '\t' '\n' | python scripts/make_sexgendertable.py > test.sexgenotype`
* --contigAnnotationsTable  
Germline contig ploidy annotations for all sex genotypes. The format is described in ContigGermlinePloidyAnnotationTableReader. For an example file, see the GATK Resource Bundle. ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/beta/GermlineCNVCaller/grch37_contig_annotations.tsv
`cat grch37_contig_annotations.tsv | sed 's/^/chr/' > hg19_contig_annotations.tsv`
* --copyNumberTransitionPriorTable  
Prior transition probabilities between different copy number states. The format is described in  IntegerCopyNumberTransitionMatrixCollection.IntegerCopyNumberTransitionMatrixCollectionReader. Lists per contig and against the sex states, the relative file paths to files that in turn list in tsv format the copy number transition priors. For an example set of files, see the GATK Resource Bundle.
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/beta/GermlineCNVCaller/copyNumberTransitionPriors/CN_transition_matrix_XX_X.tsv
`cat grch37_germline_CN_priors.tsv | sed s/^/chr/ > hg19_germline_CN_priors.tsv`
* --outputPath  
Output path for inferred model parameters, posteriors, and checkpoints.


bedtools merge -d 5 -i <(tail -n +2 targets.tsv) -c 4 -o distinct > processed_targets.bed

### CalculateTargetCoverage
https://software.broadinstitute.org/gatk/gatkdocs/4.beta.5/org_broadinstitute_hellbender_tools_exome_CalculateTargetCoverage.php
```bash
CalculateTargetCoverage --input <bamfiles> --output <outputfilename>

gatk-launch CalculateTargetCoverage --input bams/*.bam --output targetcoverage

gatk-launch CalculateTargetCoverage \
 --input bams/map-rdsbwa-D10012_8_66000_hg19s.bam \
 --input bams/map-rdsbwa-D10406_3_66001_hg19s.bam \
 --output test_targetcoverage \
 -L processed_targets.bed \
 --groupBy SAMPLE
```


### GermlineCNVCaller

```bash
 gatk-launch --javaOptions "-Xmx16g" GermlineCNVCaller \
   --jobType LEARN_AND_CALL \
   --input combined_read_counts.tsv \
   --contigAnnotationsTable grch37_contig_annotations.tsv \
   --copyNumberTransitionPriorTable grch37_germline_CN_priors.tsv \
   --outputPath learn_and_call_results \
   --sexGenotypeTable SEX_GENOTYPES.tsv \
   --disableSpark true

gatk-launch --javaOptions "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16g" GermlineCNVCaller \
 --jobType LEARN_AND_CALL \
 --input test_targetcoverage \
 --contigAnnotationsTable hg19_contig_annotations.tsv \
 --copyNumberTransitionPriorTable hg19_germline_CN_priors.tsv \
 --outputPath test_results \
 --sexGenotypeTable test.sexgenotype \
 --targets processed_targets.tsv \
 --disableSpark true \
 --rddCheckpointing false
```
