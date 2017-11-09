# GermlineCNVCaller
## Get suitable samples

## Arguments to GermlineCNVCaller
* --jobType  
either LEARN_AND_CALL or CALL_ONLY
* --input  
Combined raw or GC corrected (but not proportional) read counts table. Can be for a cohort or for a single sample. The format is described in ReadCountCollectionUtils  
* --sexGenotypeTable  
Sample sex genotypes table. The format is described in SexGenotypeTableReader, or a table listing the known sex genotypes of each sample. Only two columns, labeled SAMPLE_NAME and SEX_GENOTYPE, are required. The SEX_GENOTYPE column values must match those of the contigAnnotationsTable.
* --contigAnnotationsTable  
Germline contig ploidy annotations for all sex genotypes. The format is described in ContigGermlinePloidyAnnotationTableReader. For an example file, see the GATK Resource Bundle.
* --copyNumberTransitionPriorTable  
Prior transition probabilities between different copy number states. The format is described in  IntegerCopyNumberTransitionMatrixCollection.IntegerCopyNumberTransitionMatrixCollectionReader. Lists per contig and against the sex states, the relative file paths to files that in turn list in tsv format the copy number transition priors. For an example set of files, see the GATK Resource Bundle.
* --outputPath  
Output path for inferred model parameters, posteriors, and checkpoints.


### CalculateTargetCoverage
```bash
CalculateTargetCoverage --input <bamfiles> --output <outputfilename>
```

### ReadCountCollectionUtils

### SexGenotypeTableReader

### ContigGermlinePloidyAnnotationTableReader

### IntegerCopyNumberTransitionMatrixCollection.IntegerCopyNumberTransitionMatrixCollectionReader

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
```
