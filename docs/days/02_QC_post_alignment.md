# Quality Control Post-Alignment

## Overview
After filtering our BAM files, we need to assess the quality of our ATAC-seq experiment using specialized metrics.

### Learning Objectives
- Understand key ATAC-seq quality metrics
- Use ATAQV to generate comprehensive QC reports
- Interpret QC results to assess experiment quality



In order to assess the quality of the ATAC-seq experiment, there are several metrics that can be calculated from the aligned reads.

### 1. Key ATAC-seq Quality Metrics

**Fragment size distribution**: The distribution of fragment sizes can indicate the quality of the library preparation. A good ATAC-seq library should have a characteristic pattern with peaks corresponding to nucleosome-free regions (NFRs) and mono-, di-, and tri-nucleosomes.

**TSS enrichment**: The enrichment of reads at transcription start sites (TSS) is a key indicator of data quality. High-quality ATAC-seq data should show a strong enrichment of reads at TSSs.

### 2. ATAQV Tool

[`ATAQV`](https://github.com/ParkerLab/ataqv) is a specialized tool that calculates a variety of QC metrics for ATAC-seq data. It provides a comprehensive report that includes fragment size distribution, TSS enrichment, and other important metrics.

**Task 1: Generate individual QC reports with ATAQV**

- Create a new folder for QC results: `results/02_QC_post_alignment`  
- Run ATAQV on each filtered BAM file to generate individual QC metrics   
- Each sample will produce its own detailed QC report 


!!! ATAQV
    You can find further information about ATAQC tool and commands in the Usage section [here](https://github.com/ParkerLab/ataqv).  
    Ataqv needs a Transcription Start Site (TSS) reference file to compute TSS enrichment score. You will find this file in: `/data/references/ENCODE_mm10_M21_TSS_reference.bed`.  
    **Source**: TSS reference downloaded from [ENCODE Project](https://www.encodeproject.org/files/ENCFF498BEJ/)


??? info "Hint"
    Run `ataqv` for each `*qc_bl_filt.sorted.bam` file
    
    Key parameters needed:
    - `--tss-file`: Path to TSS reference file
    - `--metrics-file`: Output JSON file path
    - `--name`: Sample name
    - `mouse`: Genome reference




??? success "Solution"
    ```{bash}
    mkdir -p results/02_QC_post_alignment
    TSS_bed="/data/references/ENCODE_mm10_M21_TSS_reference.bed"

    for bam in results/01_filtered_bams/*qc_bl_filt.sorted.bam; do
        echo "Processing file: $bam"
        sample_name=$(basename "$bam" .qc_bl_filt.sorted.bam) 
        ataqv --name $sample_name --metrics-file results/02_QC_post_alignment/$sample_name.ataqv.json --tss-file $TSS_bed mouse $bam > results/02_QC_post_alignment/$sample_name.ataqv.out
    done
    ```
    **Parameter explanations**:
    - `--name`: Sample name for output files  
    - `--metrics-file`: Output file for metrics in JSON format  
    - `--tss-file`: BED file with TSS locations  
    - `mouse`: Genome reference (mouse or human)  
    - `$bam`: Input BAM file


**Task 2: Create multi-sample summary report**

Compile all individual QC reports into a comprehensive summary. Use the `mkarv` tool from ATAQV to combine all JSON files into an interactive HTML report

??? info "Hint"
    Use `mkarv` command to process all `*.json` files
    
    - Specify output directory name for the HTML report
    - Include all JSON files with wildcard pattern

??? success "Solution"
    ```bash
    # Create summary report from all JSON files without changing directories
    mkarv results/02_QC_post_alignment/summary_ataqv results/02_QC_post_alignment/*.ataqv.json
    ```

    **What this does**:  
    - `mkarv`: Tool to create interactive HTML summary report.  
    - `results/02_QC_post_alignment/summary_ataqv`: Output directory for the HTML report.  
    - `results/02_QC_post_alignment/*.ataqv.json`: Include all JSON files using full path.  

### 3. Interpret QC Results

Open the QC report located at: `results/02_QC_post_alignment/summary_ataqv/index.html`


**Task 3: Assess the quality of ATAC-seq experiment**

1. Overall Assessment: Do you think the experiment worked well? 
2. Fragment Size Distribution: Do you see the expected nucleosome pattern?
3. TSS Enrichment: Are the TSS enrichment scores acceptable?
4. Concerning Metrics: Are there any metrics that would concern you?

