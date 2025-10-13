# Quantify Peak Accessibility


## Overview
After building a consensus peak set, we need to quantify the number of reads overlapping each peak in each sample. This creates a count matrix for downstream differential accessibility analysis.

### Learning Objectives
- Convert peak coordinates to appropriate format for read counting
- Use featureCounts to quantify reads in peaks across samples
- Generate comprehensive quality control reports


## 1. Count reads on peaks

After defining a consensus peak set, we will count the number of reads that overlap each peak in each sample. This step can be seen as counting reads in annotated genes for RNA-seq data, where instead of genes, we are counting reads in peaks.

There are different tools that can be used for this purpose, such as [featureCounts](subread.sourceforge.net/featureCounts.html) from the Subread package or [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html).  

In this example, we will use `featureCounts`, which is a widely used tool for counting reads in genomic features.  

### 1.1 Prepare Peak Annotation for Quantification


FeatureCounts requires annotation in GTF or SAF (Simplified Annotation Format).  
The SAF format is a tab-delimited file with the following columns:   
`GeneID, Chr, Start, End, Strand`

However, your consensus peak is in BED format with only 3 columns:  

```bash
# Check file format
head results/04_consensus_peaks/consensus_peaks.bed
```

**Task1: Compare BED to SAF format**

- Which are the 2 columns you are missing to convert your BED file into SAF format? 
- Can you think of a way to provide values to those fields?

??? Answer
    `GeneID` and `Strand` fields are missing.  
    Strand is missing since the dataset has no strand information.  
    GeneID would correspond to Interval_ID, but this information has been lost while creating the consensus set, given that peaks have been merged.  
    
    `Strand`: Information is not mandatory, you can provide a `.` instead.  
    `GeneID`: You can generate sequential Interval_ID with incremental numbers (e.g., "Interval_1", "Interval_2").  



**Task 2: Convert Consensus Peaks to SAF Format**

- Create a new folder named: `results/05_counts_reads`
- Use awk or an alternative way to convert the bed file to SAF format. Don't forget to add a header "GeneID    Chr Start   End Strand". 

??? info "Hint"
    Use `awk` to:  
    - Add sequential peak IDs (e.g., "Interval_1", "Interval_2").  
    - Convert BED coordinates to SAF format.  
    - Add strand information (use "." for unstranded).  

??? success "Solution"

    ```{bash}
    # Create output directory
    mkdir results/05_counts_reads

    # Create SAF file with header
    echo "GeneID    Chr Start   End Strand" > results/04_consensus_peaks/consensus_peaks.saf

    # Convert BED to SAF format
    awk '{OFS = "\t"} {print "Interval_"NR,$1,$2,$3,"."}' results/04_consensus_peaks/consensus_peaks.bed >> results/04_consensus_peaks/consensus_peaks.saf

    # check the new file format
    head results/04_consensus_peaks/consensus_peaks.saf
    ```
    **What this does:**.  
    - `echo -e`: Creates header with tab separators.  
    - `awk OFS="\t"`: Sets output field separator to tab.  
    - `"Interval_"NR`: Creates unique ID using row number.  
    - `$1,$2,$3`: BED coordinates (chr, start, end).  
    - `"."`: Unstranded (appropriate for ATAC-seq peaks).  

Now we can use featureCounts to count the reads in each peak for each sample

**Task 3: Count reads in peaks**

- Run featureCounts to quantify how many filtered reads (in `*qc_bl_filt.sorted.bam`) overlap consensus peaks in each sample.  
- Use paired-end counting mode
- Count all filtered BAM files simultaneously

!!! FeatureCounts
    `-F SAF`: specify that the annotation file is in SAF format  
    `-p`: specify that the input files are paired-end  
    `-a`: specify the annotation file  
    `-o`: specify the output file  
    
    These are only some of all parameteres one can apply. You can also adjust the parameters of featureCounts based on your specific requirements, such as setting a minimum mapping quality or handling multi-mapping reads.


??? success "Solution"

    ```{bash}
    path_bams="results/01_filtered_bams"
    featureCounts -F SAF -T 2 -p -a results/04_consensus_peaks/consensus_peaks.saf -o results/05_counts_reads/feature_Counts.txt $path_bams/*qc_bl_filt.sorted.bam 2> results/05_counts_reads/featureCounts.log
    ```

The output file will contain the counts of reads in each peak for each sample, which can be used for downstream analysis such as differential accessibility analysis.

**Task 4: Examine the results**

Have a look to the output, which file contains the counts matrix?:
```{bash}
ls results/05_counts_reads/

```

??? Answer
    FeatureCounts generates:  
    - `feature_Counts.txt`: Main count matrix.  
    - `feature_Counts.txt.summary`: Counting statistics.  
    - `featureCounts.log`: Processing log and any warnings.  
  
    File `results/05_counts_reads/feature_Counts.txt` contains the counts table with:  
    - Rows: Annotated peaks.  
    - Columns: peaks annotation information (from SAF file) + Samples (counts in each sample)


- Check how many peaks contains the counts matrix, is it the same number you had in your consensus peak set?

```bash
# Count number of peaks quantified
wc -l results/05_counts_reads/feature_Counts.txt
```
??? Answer
    Yes, it is the same number of peaks, although not same number of rows because `feature_Counts.txt` has 2 lines header.


## 2.  Quality Control Assessment

After filtering, peak calling and counting reads in peaks, we can run MultiQC to aggregate the QC metrics from different steps and generate a comprehensive report.  

This will help us to assess further the overall quality of the ATAC-seq data and identify any potential issues that may need to be addressed before proceeding with downstream analysis.  

**Task 3: Generate MultiQC Report**  

- Run multiqc on the entire `results` directory


??? success "Solution"

    ```{bash}
    multiqc --outdir results/multiQC_report --title "ATAC-seq_Pipeline_Summary" results/

    ```
    **Parameter explanations:**.  
    - `--outdir`: Output directory for the report.  
    - `--title`: Custom title for the report.  
    - `results/`: Search this directory for QC files.  


Download and open the report and explore the QC metrics:

- How many peaks were successfully quantified?   
- What percentage of reads were assigned to peaks in each sample?   
- Are there differences in counting efficiency between samples?   
- Do any samples show concerning quality metrics?   