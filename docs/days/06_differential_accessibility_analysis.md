# Differential Accessibility Analysis



## Overview
After quantifying reads in consensus peaks, we will identify regions with significantly different accessibility between conditions using DESeq2.

### Learning Objectives
- Import and prepare count data for differential analysis
- Create proper metadata for experimental design
- Perform differential accessibility analysis with DESeq2
- Visualize and interpret results
- Prepare data for downstream functional analysis



!!! Note 
    For this section, you have to move to **R**. Change from `Terminal` to `Console` in your Rstudio server. To do that, click on the **`Console`** tab on the top left of your terminal. 

Load necessary packages for the analysis
    ```r
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    library(GenomicRanges)
    ```

## 1. Import read counts and metadata


Start by loading featureCounts output as a matrix of counts


**Task 1: Load counts matrix**

- Read the output of featurecounts that contains counts for all peaks in all samples: `results/05_counts_reads/feature_Counts.txt`
- Keep the header of the table
- Process the table in order to keep: read-counts columns only, peak IDs as row.names, simple column names (sample names)
- Convert the table into matrix format

??? success "Solution"

    ```r
    # Read table of counts (output from FeatureCounts)
    counts_file <- "results/05_counts_reads/feature_Counts.txt"
    cts_table <- read.table(counts_file, header = T)
    
    # convert into matrix, with interval IDs as rownames
    cts <- cts_table[,7:10]
    row.names(cts) <- cts_table[,1]
    colnames(cts) <- gsub(colnames(cts), pattern = "results.01_filtered_bams.", replacement = "") %>% 
                     gsub(., pattern=".qc_bl_filt.sorted.bam", replacement="")
    #check the object 
    cts
    ```

**Task 2: Create sample metadata**

- Create metadata table as data.frame
- Each row is a sample, columns contain metadata information about samples (e.g. condition/treatment/batch...). Sample names are row.names
- Here you only have information about condition: Kidney and Cerebrum. The metadata will contain only one column. If you would have a second factor, e.g. treatment (treated vs untreated), you would need to add this information as a second column in the data.frame.  
- Use factor variables (not character variables)

!!! Warning
  
    Important: The order of the rows in the data.frame must match the order of the columns in the counts matrix
               The metadata data.frame must have factor variables (not character variables)


```r
condition <- factor( c(rep("Cerebrum",2), rep("Kidney",2)) )
colData <- data.frame(condition, row.names = colnames(cts))

# Verify the setup
print(colData)
all(rownames(colData) == colnames(cts))
```

## 2. DESeq2 Differential Analysis

Bring together counts and metadata to create a DESeq object

```r
dds <- DESeqDataSetFromMatrix(
  countData = cts, colData = colData, 
  design = ~ condition)
dim(dds)
```

Optional: (Remove peaks with insufficient read counts)

```r
idx <- rowSums(counts(dds, normalized=FALSE) >= 30) >= 2
dds.f <- dds[idx, ]
dim(dds.f)
```

Perform the estimation of dispersions

```r
dds <- DESeq(dds)
```

And plot PCA of the samples
```r
# we apply the variance stabilising transformation (vst) to make the read counts comparable across libraries
# (nb : this is not needed for DESeq DE analysis, but rather for visualisations that compare expression across samples, such as PCA. This replaces normal PCA scaling)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("condition"), ntop="all")
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))
plotPCA(vsd, intgroup=c("sizeFactor"), ntop="all")
```


## 3. Explore DA results. 

After the dispersion estimates have been calculated, you can proceed to test for differential accessibile (DA) regions between the two conditions (Kidney vs Cerebrum).  

The `results()` function extracts the results table from DA analysis with the log2 fold changes, p-values and adjusted p-values for each peak.  

```r
DA_results <- results(dds)
summary( DA_results )
head(DA_results)
```

Have a look in `DA_results` object.  

- How many peaks have been analysed?
- At a cutoff of padj < 0.01, how many significant DA peaks there are?
- How many have singificantly higher accessibility in Kidney compared to Cerebrum? 


??? Success "Solution" 
    ```r
        cat("Total peaks analyzed:", nrow(DA_results), "\n")
        cat("Significant peaks:", sum(DA_results$padj < 0.01, na.rm = TRUE), "\n")
        cat("Upregulated in Kidney:", sum(DA_results$log2FoldChange > 0 & 
                                        DA_results$padj < 0.01, na.rm = TRUE), "\n")
        cat("Upregulated in Cerebrum:", sum(DA_results$log2FoldChange < 0 & 
                                        DA_results$padj < 0.01, na.rm = TRUE), "\n")
    ```

You can visualise the results of DA peaks with a Volcano plot, analogous to RNAseq DE genes results

```r
 # Create volcano plot
    volcano_data <- as.data.frame(DA_results)
    volcano_data$significant <- abs(volcano_data$log2FoldChange) > 2 & volcano_data$padj < 0.01
    
    ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = significant), alpha = 0.6) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", alpha = 0.5) +
      labs(x = "Log2 Fold Change", 
           y = "-Log10 Adjusted P-value",
           title = "Volcano Plot: Kidney vs Cerebrum",
           subtitle = "Red points: |Log2FC| > 2 and padj < 0.01") +
      theme_minimal()
```

!!! Note

    By default, the results() function will extract the results for the last variable in the design formula (here: condition) and will perform a comparison of the second level of the factor over the first level (here: Kidney over Cerebrum).
    If you want to extract results for a different comparison, you can specify the contrast argument.  
    Have a look at [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) documentation for more information.



## 4. Save differential accessibility results for downstream analysis:

Save differential accessibility results for downstream analysis:


```r
dir.create("results/06_DA_analysis")
write.table(DA_results, file="results/06_DA_analysis/DA_results.txt", quote=FALSE)
```




## 5. Create GRanges object


Convert the peak coordinates into a GRanges object (from [GenomicRanges package](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html)).  
`GRanges` class represents a collection of genomic ranges and associated data to them. 
In this case, you will use it to represent the peak coordinates, and add as metadata the log2 fold changes and adjusted p-values from the differential accessibility analysis.


```r
    # Prepare peak coordinates object 
    peaks_coord <- cts_table[,1:4]
    head(peaks_coord)

    # Select information from DESeq2 you want to keep
    DA_stats <- as.data.frame(DA_results)[c(2,3,5,6)]
    head(DA_stats)

    # Merge both objects and convert them into GRanges class
    peaks_df <- merge(DA_stats, peaks_coord, by.x="row.names", by.y="Geneid")
    peaks_df$Chr <- paste0("chr",peaks_df$Chr)
    peaks_gr = makeGRangesFromDataFrame(peaks_df, keep.extra.columns=T)

    # Have a look into the new object
    head(peaks_gr)
```
To this object `peaks_gr`, you will add other metadata along the downstream analysis, like genomic overlap.

Save the object
```r
saveRDS(peaks_gr, "results/06_DA_analysis/DA_results.RDS")
```

And clean the environment

!!! Warning
    You must save the objects before running the following code. After running it your variables and objects will be removed.

```r
rm(list = ls())
.rs.restartR()
```