# GermlineCNVCaller
## Get suitable samples

## Arguments to GermlineCNVCaller
* --jobType  
either LEARN_AND_CALL or CALL_ONLY
* --input  
Combined raw or GC corrected (but not proportional) read counts table. The format is described in ReadCountCollectionUtils  
Can be for a cohort or for a single sample.
* --sexGenotypeTable  
Sample sex genotypes table. The format is described in SexGenotypeTableReader
* Germline contig ploidy annotations for all sex genotypes (specified via argument --contigAnnotationsTable). The format is described in ContigGermlinePloidyAnnotationTableReader
* Prior transition probabilities between different copy number states (specified via argument --copyNumberTransitionPriorTable). The format is described in  IntegerCopyNumberTransitionMatrixCollection.IntegerCopyNumberTransitionMatrixCollectionReader
* Output path for inferred model parameters, posteriors, and checkpoints (specified via argument --outputPath)

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
