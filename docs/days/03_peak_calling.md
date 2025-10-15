# Peak Calling

## Overview
After quality control, we will identify accessible chromatin regions (peaks) using different approaches and fragment size filters.  

We will compare peak calling results using:  
1. Only nucleosome-free fragments  
2. All fragments combined  
3. HMMRATAC (ATAC-seq specific peak caller)  

### Learning Objectives
- Understand nucleosome-free vs. nucleosome-associated fragments
- Use MACS3 for peak calling with different fragment size filters
- Compare peak calling results from different approaches
- Learn about ATAC-seq specific peak caller HMMRATAC  


## 1. Fragment size filtering

Before calling peaks with MACS3, we will separate fragments based on their size:

- **Nucleosome-free (NF) fragments** (< 100 bp): Represent open chromatin regions
- **Nucleosome-associated fragments** (> 100 bp): Represent regions with positioned nucleosomes


**Task 1: Filter for nucleosome-free (NF) fragments**

We can do that using [`samtools view`](https://www.htslib.org/doc/samtools-view.html) as we have seen in previous exercises, and filter reads based on fragment length.  


- Which column from the .bam file contains fragment length information?  

??? info "Answer"
    In SAM format, column 9 is described as "TLEN: signed observed Template LENgth", which corresponds to the insert size length. 

With the following code, you will:  
- Create output directory: `results/01_filtered_bams/NF_bams`.  
- Use `samtool view` to filter BAM files to keep only fragments with length 1-100 bp.    
- Save filtered files as: `results/01_filtered_bams/NF_bams/${sample_name}_NF.bam`.  
- Maintain BAM file headers during filtering.  


!!! info "Code"

    ```bash
    # create new directory for peaks
    mkdir -p results/01_filtered_bams/NF_bams

    # filter for NF fragments only

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        # run command
        samtools view -h $path_bams/$sample_name.qc_bl_filt.sorted.bam | awk 'substr($0,1,1)=="@" || ($9>= 1 && $9<=100) || ($9<=-1 && $9>=-100)' | \
        samtools view -b > $path_bams/NF_bams/${sample_name}_NF.bam 
    done 
    ```

    **Parameter explanations**:  
    - `-h`: Keep header in output  
    - `-b`: Output in BAM format  
    - `$9`: Insert size (TLEN field in SAM format)  
    - **Filter logic**: Keep fragments with insert size between 1 and 100bp (nucleosome-free regions)


**Task 2: sort and index BAM files**

- Sort and index the bams files for next step
- Remove the unsorted file


```bash
for sample_name in "${samples[@]}"; do
    echo "Processing sample: $sample_name"

    # run command
    samtools sort -o $path_bams/NF_bams/${sample_name}_NF.sorted.bam $path_bams/NF_bams/${sample_name}_NF.bam
    samtools index $path_bams/NF_bams/${sample_name}_NF.sorted.bam
done
```

!!! Warning
    **Important before deleting:**  
    Run this command *only* if you have named your files exactly as specified in this tutorial and have completed all sorting and indexing steps from the previous command.

Remove intermediate files
```bash
rm results/01_filtered_bams/NF_bams/*_NF.bam
```

## 2. Peak calling strategies


### 2.1 MACS3 Peak Calling

MACS3 is a widely used peak calling tool that can handle both narrow and broad peaks. We'll start using MACS3 function "callpeak" and compare results using different fragment size filters.

#### Peak Calling with Nucleosome-Free Fragments

We will first call peaks using MACS3 on NF reads only, focusing on fragments most likely to represent open chromatin regions.


**Task 3: MACS3 on NF fragments**

- Create a new folder named: `results/03_peak_calling`
- Create a subfolder inside called: `NF_peaks`
- Call peaks on NF reads using MACS3 callpeak function. Save the results inside `results/03_peak_calling/NF_peaks/` and name the files as: `NF_peaks_${sample_name}`


!!! MACS3 parameters
    You can have a look at the MACS3 documentation for more details on the parameters used [here](https://macs3-project.github.io/MACS/docs/callpeak.html)

??? info "Hint"
    For ATAC-seq data, we will use the BAMPE format, which is suitable for paired-end data and we will set the genome size to "mm" for mouse. We will also set a q-value cutoff of 0.01 to control the false discovery rate.


??? success "Solution"

    ```bash
    mkdir -p results/03_peak_calling/NF_peaks
    path_peaks="results/03_peak_calling/"

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        # run command
        macs3 callpeak -f BAMPE -t $path_bams/NF_bams/${sample_name}_NF.sorted.bam -g mm -q 0.01 --name ${sample_name}_NF --outdir results/03_peak_calling/NF_peaks/NF_peaks_${sample_name}/
    done
    ```
    parameters explanation:  
    `-f BAMPE`: input file format is BAM paired-end  
    `-t`: input file (BAM file)  
    `-g mm`: It’s the mappable genome size or effective genome size (some are pre-computed, like mouse, and you can specify "mm"
    `-q 0.01`: q-value cutoff for peak detection  
    `--outdir`: output directory for peak files  


#### Peak Calling with all Fragments

We will do the same, but using all fragments 

**Task 4: MACS3 on all fragments**

- Create a new folder named: `results/03_peak_calling/all_peaks`
- Call peaks on all filtered reads using MACS3 callpeak function. Save the results inside `results/03_peak_calling/all_peaks/` and name the files as: `all_peaks_${sample_name}`

??? success "Solution"

    ```bash
    mkdir -p results/03_peak_calling/all_peaks

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample_name"
        
        # run command
        macs3 callpeak -f BAMPE -t $path_bams/${sample_name}.qc_bl_filt.sorted.bam -g mm -q 0.01 --name ${sample_name}_all --outdir $path_peaks/all_peaks/all_peaks_${sample_name}/ 2> $path_peaks/all_peaks/${sample_name}_macs3.log
    done
    ```



### 2.2 HMMRATAC Peak Calling

[HMMRATAC](https://academic.oup.com/nar/article/47/16/e91/5519166?login=true) is specifically designed for ATAC-seq data and uses a Hidden Markov Model to identify accessible chromatin regions. It models the distribution of fragment lengths to distinguish between nucleosome-free regions and nucleosome-bound regions, to provide a more accurate identification of open chromatin regions.  

We will use the filtered BAM files (all fragments) for this analysis.  

!!! tip "Bonus Task"

    - Create a new directory named: `results/03_peak_calling/hmmratac_peaks`.
    - Have a look at the [HMMRATAC documentation](https://macs3-project.github.io/MACS/docs/hmmratac.html) for more details on the parameters used.
    - This step can take a bit longer, while waiting you can move to the next section of the tutorial.

    ??? success "Solution"

        ```bash
        mkdir -p results/03_peak_calling/hmmratac_peaks

        for sample_name in "${samples[@]}"; do
            echo "Processing sample: $sample_name"
            macs3 hmmratac -i $path_bams/${sample_name}.qc_bl_filt.sorted.bam -f BAMPE --name ${sample_name}_hmmratac --outdir $path_peaks/hmmratac_peaks/hmmratac_peaks_${sample_name}/
        done
        ```
        Parameters explanation:  
        `-i`: input file (BAM file)  
        `-f BAMPE`: input file format is BAM paired-end  
        `--outdir`: output directory for peak files  


## 3. Compare Peak Calling Results

Now let's visualise and compare the different peak calling results using the Integrative Genomics Viewer (IGV) to understand how each method performs.

**Task 5: Load traks into IGV** 

For **Cerebrum_rep1** sample, download and load the following files in IGV:

#### Peak Files:
- MACS3 NF peaks: `results/03_peak_calling/NF_peaks/NF_peaks_Cerebrum_rep1/Cerebrum_rep1_NF_peaks.narrowPeak`
- MACS3 all peaks: `results/03_peak_calling/all_peaks/all_peaks_Cerebrum_rep1/Cerebrum_rep1_all_peaks.narrowPeak`
- HMMRATAC peaks: `results/03_peak_calling/hmmratac_peaks/hmmratac_peaks_Cerebrum_rep1/Cerebrum_rep1_hmmratac_accessible_regions.narrowPeak` *(if you completed the bonus task, if you didn't and you are curious, you will find the files in `/data/Solutions/`)*

#### BAM Files for Read Coverage:
- All fragments BAM: `results/01_filtered_bams/Cerebrum_rep1.qc_bl_filt.sorted.bam`
- Nucleosome-free fragments BAM: `results/01_filtered_bams/NF_bams/Cerebrum_rep1_NF.sorted.bam`

**Task 6: Compare results**  

- Do you observe big differences between methods?  
- Which method would work better to study TF binding sites? and for DA analysis?   
- Can you find an example of a gene where the nucleosomal pattern at TSS can be percieved?  
