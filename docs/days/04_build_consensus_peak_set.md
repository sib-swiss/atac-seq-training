# Build Consensus Peak Set

## Overview
After calling peaks in individual replicates, we need to create a unified peak set for downstream analysis by identifying reproducible peaks across replicates.

### Learning Objectives
- Understand the importance of reproducible peaks in ATAC-seq analysis
- Use bedtools to intersect peaks between replicates
- Create a consensus peak set for differential accessibility analysis
- Visualize consensus peaks in IGV



## 1. Build Consensus Peak Annotation

For robust downstream analysis, we need peaks that are reproducible across biological replicates. This ensures our analysis focuses on genuine accessible chromatin regions rather than technical artifacts.

Build a consensus peak set, which will become your reference annotation for counting reads in peaks later on.


**Task 1: intersect peaks between replicates**

Find peaks that are present in both replicates within each condition:

- Create output directory: `results/04_consensus_peaks`
- Use bedtools intersect with 25% reciprocal overlap requirement
- Generate separate intersection files for each condition
- In order to get original peaks from both replicates (without appending the output), run bedtools intersect twice, swapping the order to capture all peaks. 

??? info "Hint"
    Use `bedtools intersect` with parameters:  
    - `-f 0.25`: minimum 25% overlap  
    - `-r`: reciprocal overlap required  
    - `-wa`: write original entries from file A  
    - Run twice with swapped files to capture all overlapping peaks

??? success "Solution"
    ```bash
    # Create output directory
    mkdir -p results/04_consensus_peaks
    path_peaks="results/03_peak_calling/all_peaks"

    # Intersect Kidney replicates (run twice to get all overlaps)
    bedtools intersect -wa -a $path_peaks/all_peaks_Kidney_rep1/Kidney_rep1_all_peaks.narrowPeak \
                          -b $path_peaks/all_peaks_Kidney_rep2/Kidney_rep2_all_peaks.narrowPeak \
                          -f 0.25 -r > results/04_consensus_peaks/Kidney_intersect.bed
    
    bedtools intersect -wa -b $path_peaks/all_peaks_Kidney_rep1/Kidney_rep1_all_peaks.narrowPeak \
                          -a $path_peaks/all_peaks_Kidney_rep2/Kidney_rep2_all_peaks.narrowPeak \
                          -f 0.25 -r >> results/04_consensus_peaks/Kidney_intersect.bed

    # Intersect Cerebrum replicates (run twice to get all overlaps)
    bedtools intersect -wa -a $path_peaks/all_peaks_Cerebrum_rep1/Cerebrum_rep1_all_peaks.narrowPeak \
                          -b $path_peaks/all_peaks_Cerebrum_rep2/Cerebrum_rep2_all_peaks.narrowPeak \
                          -f 0.25 -r > results/04_consensus_peaks/Cerebrum_intersect.bed
    
    bedtools intersect -wa -b $path_peaks/all_peaks_Cerebrum_rep1/Cerebrum_rep1_all_peaks.narrowPeak \
                          -a $path_peaks/all_peaks_Cerebrum_rep2/Cerebrum_rep2_all_peaks.narrowPeak \
                          -f 0.25 -r >> results/04_consensus_peaks/Cerebrum_intersect.bed

    # Check file format
    head results/04_consensus_peaks/Cerebrum_intersect.bed
    ```

    **Parameter explanations**:  
    - `-a`: First input file (bed or narrowPeak format)  
    - `-b`: Second input file (bed or narrowPeak format)    
    - `-f 0.25`: Minimum overlap required as a fraction of A  
    - `-r`: Require reciprocal overlap  
    - `-wa`: Write the original entry in A for each overlap  
    - **Why twice?**: Running bedtools intersect twice with swapped files ensures we capture all overlapping peaks  


**Task 2: Create final consensus peak set**

Merge the intersected peaks from both conditions into a final consensus peak set, allowing a maximum distance of 10bp between peaks to be merged (-d 10).  

- Combine peaks from both conditions
- Sort by genomic coordinates
- Merge nearby peaks (within 10bp) to avoid redundancy




```{bash}
# concatenate intersection of both replicates and sort the file
cat results/04_consensus_peaks/Kidney_intersect.bed results/04_consensus_peaks/Cerebrum_intersect.bed | sort -k1,1 -k2,2n > results/04_consensus_peaks/consensus_peaks_temp.bed
sort -k1,1 -k2,2n results/04_consensus_peaks/consensus_peaks_temp.bed > results/04_consensus_peaks/consensus_peaks_temp_sorted.bed

# Use bedtools merge to merfe peaks at shorter distance than 10bp
bedtools merge -d 10 -i results/04_consensus_peaks/consensus_peaks_temp.bed > results/04_consensus_peaks/consensus_peaks.bed
```


Add the `results/04_consensus_peaks/consensus_peaks.bed` track to IGV and have a look to the resulting peak annotation. 

**Task 3: Visualize Consensus Peaks**


- Load `results/04_consensus_peaks/consensus_peaks.bed` in IGV
- Compare with individual replicate peak files (for that, load results from peak calling from each sample)
- Assess the quality and coverage of your consensus peak set

1. Do the consensus peaks cover the main accessible regions you observed in individual replicates?
2. Was the merging of nearby regions to stringent? too relax? 
3. Are there regions where you lost peaks due to the stringent overlap requirements?


**Task 4: Count number of peaks**

- Count how many peaks contains your consensus peak set.
- Count how many peaks each individual sample contained

```bash
# Count final consensus peaks
wc -l results/04_consensus_peaks/consensus_peaks.bed

# Count peaks in each sample
samples=(Kidney_rep1 Kidney_rep2 Cerebrum_rep1 Cerebrum_rep2)
echo "Number of peaks called by MACS3:"
for sample in "${samples[@]}"; do
    wc -l results/03_peak_calling/all_peaks/all_peaks_${sample}/NA_peaks.narrowPeak
done

# Check file format
head results/04_consensus_peaks/consensus_peaks.bed
```