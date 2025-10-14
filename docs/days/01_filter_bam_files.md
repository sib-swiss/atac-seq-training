# Post-Alignment Reads Filtering

## Overview
In this section, we will filter BAM files to remove low-quality reads and prepare them for downstream analysis.

### Learning Objectives
- Apply quality filters to BAM files using samtools
- Remove reads from blacklisted genomic regions
- Understand ATAC-seq specific filtering considerations


## 0. Dataset

We are going to work with a subset of the publicly available ATAC-seq dataset from [Liu et al. 2019](https://www.nature.com/articles/s41597-019-0071-0). We will process and compare data from 2 different adult mouse tissues: 

- **Kidney**: Rep1, Rep2     
- **Cerebrum**: Rep1, Rep2  

  
In order to save time, raw reads have already been trimmed using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and mapped to the reference genome (mm10) using [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) in end-to-end mode. Therefore, we will start the analysis directly from the alignment output (**.bam** files).  

To avoid large waiting times during high-demanding computational steps, .bam files have been subset to keep only reads aligning to chromosome 6.  

!!! Note

    You can find the .bam files in `/data/Liu_alignments_chr6/` folder.  
    Raw data can be found in [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/SRP167062).



## 1. Filtering Reads from Alignments

### 1.1 Common NGS Filtering Steps

To ensure high-quality data, it's important to filter aligned reads based on standard NGS criteria, such as:

- Removing duplicate reads (often PCR artifacts).  
- Applying minimum mapping quality thresholds.  

These filters help reduce noise and improve the accuracy of downstream analyses.

!!! Note
    Tools like Picard can identify and mark duplicate reads. In this dataset, duplicates have already been marked in the provided .bam files. To learn more about Picard have a look [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).   
   

A useful tool for filtering alignments is [`samtools view`](https://www.htslib.org/doc/samtools-view.html). It allows printing all the alignments from a .bam file in SAM format, with the option to previously filter them based on specific flags (e.g. mapped/unmapped, primary/secondary, duplicates).  

In the next exercises, you will use samtools view to visualise the alignment output and filter .bam files according to various criteria.


**Task 1: Visualise alignment output** 

- Visualise one of the .bam files with samtools running the following command:

```bash
samtools view /data/Liu_alignments_chr6/Kidney_rep1.bam | head
```
Can you recognise what the different columns represent?  
Why the 3rd column contains only number 6?  
What is the 2nd column? 

??? done "Answer"
    The columns correspond to SAM format. You can find an explanation of each column [here](https://www.htslib.org/doc/sam.html).  
    The 3rd column represents the chromosome name. In this tutorial, .bam files have been subset to contain only chr6, thus it makes sense they are all the same.  
    The 2nd column contains the FLAG, you can use this value to filter out or in certain reads.  


**Task 2: Quality filtering with samtools** 

Check in this [website](https://broadinstitute.github.io/picard/explain-flags.html), which FLAG number  would you need in order to filter out from the alignments:  

- Unmapped reads
- Reads who's mate is unmapped
- Read duplicates
- Not primary alignment reads

??? done "Answer"
    You need flag: 1292


Now, use samtools view [manual](https://www.htslib.org/doc/samtools-view.html) to filter the previous .bam file based on all these criteria:  

- Filter out reads with flag 1292.  
- Keep only paired reads.  
- Filter out reads with mapping quality < 10.  

<details>
<summary>Hint</summary>
You can use:   
    -f and -F flags to filter in or out reads (respectively) based on a combination of SAM Flags; use -q: for MapQ threshold
</details>

    
??? success "Solution"

    ```bash
    samtools view  -q 10 -f 1 -F 1292 /data/Liu_alignments_chr6/Kidney_rep2.bam | head
    ```

Finally, apply these filters to all .bam files and save them in .bam format (check parameter -b for that) in `results/01_filtered_bams/<sample_name>.qc_filt.bam`





??? success "Solution"

    ```bash
    # To be able to iterate through all samples:
    samples=(Kidney_rep1 Kidney_rep2 Cerebrum_rep1 Cerebrum_rep2)

    # create new directory for filtered bams
    mkdir -p results/01_filtered_bams

    # save the path to the folder in a variable
    path_bams="results/01_filtered_bams"


    # filter files using a for loop

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        # run command
        samtools view -h -b -q 10 -f 1 -F 1292 -o $path_bams/$sample_name.qc_filt.bam /data/Liu_alignments_chr6/$sample_name.bam
    done
    ```

    - **-h**: keep header  
    - **-b**: output in bam format  
    - **-q 10**: minimum mapping quality 10  
    - **-f 1**: keep paired  
    - **-F 1292**: exclude reads with any of the following flags: read unmapped, mate unmapped, not primary alignment, read is duplicate


After filtering we will sort and index the bam files for the next step


**Task 3: Sort and index BAM files** 

Next, sort and index the bam files for downstream analysis.

```bash
for sample_name in "${samples[@]}"; do
    echo "Processing sample: $sample_name"
    
    # run commands
    samtools sort -o $path_bams/$sample_name.qc_filt.sorted.bam $path_bams/$sample_name.qc_filt.bam
    samtools index $path_bams/$sample_name.qc_filt.sorted.bam
done
``` 


To avoid confusion and large size files, keep only the sorted and indexed bams:

!!! Warning
    **Important before deleting:**  
    Run this command *only* if you have named your files exactly as specified in this tutorial and have completed all sorting and indexing steps from the previous command.

```bash
rm results/01_filtered_bams/*qc_filt.bam
```


### 1.2 ATACseq related filtering steps

Mitochondria DNA is nucleosome free, therefore it is more accessible for Tn5 and several reads may have originated from mitochondrial DNA. After having assessed mitochondrial % on the QC, we can discard reads coming from chmMT to avoid biases in downstream analysis.

Since we are working only with chr6 we don't need to do this step, but here is the command you could use for that:

```bash
samtools view -h input.bam | awk  '($3 != "MT")' | samtools view -hb - > output.bam
```
!!! Note
    The mitochondrial chromosome name may differ depending on the reference genome (e.g., "MT", "chrM", "chrMT").

Next, we will remove reads overlapping problematic regions of the genome. ENCODE consortium has created comprehensive lists of such regions (anomalous, unstructured or high signal in NGS experiments) for different genome species (including mouse mm10). These lists are called ENCODE Blacklists, and you can find them [here](https://github.com/Boyle-Lab/Blacklist/). 

!!! Note
    The regions for mm10 have been dowloaded as a .bed file, you can find it here: `/data/references/mm10-blacklist.v2.nochr.bed`

**Task 4: Remove blacklist regions**


- Using [`bedtools intersect`](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) filter out reads overlapping regions in the Blacklist .bed file
- Use previously filtered .bam files as input:  
`results/01_filtered_bams/*qc_filt.sorted.bam`, don't forget to specify your input is in BAM format
- Use the `/data/references/mm10-blacklist.v2.nochr.bed` as regions to filter out reads from (it is already sorted)
- Save the results in `results/01_filtered_bams/` in BAM format with the following output name: `<sample_name>.qc_bl_filt.bam` 

!!! Bedtools intersect useful information
    `Bedtools intersect` allows one to screen for overlaps between two sets of genomic features/regions, and then decide on which kind of information do you want to report.  
    Here, we will intersect:  
    a) Aligned reads to the genome (filtered and sorted .bam files)  
    b) Problematic regions listed in Blacklist .bed file  
    We do not want to keep the reads that overlap Blacklist regions.  
    You can find documentation on which parameters to use [here](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)



??? success "Solution"  

    ```bash
    blacklist="/data/references/mm10-blacklist.v2.nochr.bed"

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        # run command
        bedtools intersect -v -abam $path_bams/$sample_name.qc_filt.sorted.bam -b $blacklist > $path_bams/$sample_name.qc_bl_filt.bam
    done
    ```

     Parameter explanation:  
     **-v**: only report those entries in A that have no overlap with B  
     **-abam**: input is in bam format 



**Task 5: Sort and index the previous files**  

Sort and index the bam files for downstream analysis.  


```bash
for sample_name in "${samples[@]}"; do
    echo "Processing sample: $sample_name"
    
    # run commands
    samtools sort -o $path_bams/$sample_name.qc_bl_filt.sorted.bam $path_bams/$sample_name.qc_bl_filt.bam
    samtools index $path_bams/$sample_name.qc_bl_filt.sorted.bam
done
```

After sorting the .bam files, we don't need the unsorted .bam. To free some space we will remove the unsorted bam files (intermediate files)
.
!!! Warning
    **Important before deleting:**  
    Run this command *only* if you have named your files exactly as specified in this tutorial and have completed all sorting and indexing steps from the previous command.

```bash

rm results/01_filtered_bams/*qc_bl_filt.bam
```


!!! tip "Bonus Task"
    **Taks 6: Convert .bam file into .bigwig format**

    In order to visualise ATAC-seq coverage, you can convert .bam files into .bigwig format, used to represent coverage tracks in genome browsers like IGV.  
    A tool that allows you to do this conversion is [`deeptools bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)

    Next, convert the sorted and filterd .bams into .bigwigs


    ```bash
    mkdir results/01_filtered_bams/filt_bigwigs
    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        #run command
        bamCoverage -b $bam -o $path_bams/filt_bigwigs/$sample_name.bw --region 6:1500000:20000000 --normalizeUsing CPM
    done
    ```


    Download the .bigwigs and load them into IGV. 
    Remember to use mm10 genome.  

    Look into Asns gene, which tissue has higher DNA accessibility? 

