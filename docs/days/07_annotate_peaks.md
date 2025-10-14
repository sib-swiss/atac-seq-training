# Peak Annotation and Functional Analysis

After performing differential accessibility analysis with DESeq2, you will annotate peaks based on their genomic location using ChIPseeker and perform functional enrichment analysis.

### Learning Objectives
- Load differential accessibility results and prepare for annotation
- Understand genomic feature annotation using ChIPseeker
- Visualize peak distribution and genomic context
- Perform functional enrichment analysis on peak-associated genes 
- Interpret biological significance of accessibility changes 


Load necessary packages for the analysis 

```r
# Load libraries
library(ChIPseeker)
require("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(clusterProfiler)
library("org.Mm.eg.db")


# Set up annotation databases
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
AnnoDb <- 'org.Mm.eg.db'
```


## 1. Load DA results

Load the DA results you have generated in the previous section, and name it `DA_results` again:

??? Success "Solution"
    ```r
    gr <- readRDS("results/06_DA_analysis/DA_results.RDS")
    head(gr)
    ```

## 2. Visualize Peak Distribution

Let's visualize how peaks are distributed across the genome. In this case, they are all located at chr6, but when analysing the entire genome you may want to see if they cluster in specific chromosomes or locations. 
`covplot` function from [ChipSeeker](https://academic.oup.com/bioinformatics/article/31/14/2382/255379) package allows you to do that: 

**Task1: Plot peaks across the genome**

- Serach for `covplot` function in [ChipSeeker documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)  and plot the genomic location of peaks, as well as their logFC. 

??? Success "Solution"
    ```{r}
    # Create coverage plot showing peak positions and their fold changes
    covplot(gr, weightCol = "log2FoldChange")
    ```

This plot shows the distribution of peaks across chromosomes, with y-axis representing the log2 fold change values.


## 3. Filter Significant Peaks

Split peaks into upregulated and downregulated categories:

- Sort GRanges by genomic coordinates 
- Filter for significant peaks (padj < 0.01)
- Apply fold change thresholds (|log2FC| > 2)

```r
# Sort the GRanges object by genomic coordinates
gr <- sort(gr, by = ~ seqnames + start + end)

# Split peaks into upregulated and downregulated based on significance and fold change
gr_list <- list(
  up = gr[gr$padj < 0.01 & gr$log2FoldChange > 2,], 
  down = gr[gr$padj < 0.01 & gr$log2FoldChange < -2,]
)

# Check the structure of gr_lists (list with 2 GRanges object)
head(gr_list)
```

**Task 2: Explore gr objects** 

- Check how many rows contain each object of the list, are the numbers balanced?

??? Success "Solution"
    ```{r}
    # Count peaks in each category
    cat("Number of upregulated peaks:", length(gr_list$up), "\n")
    cat("Number of downregulated peaks:", length(gr_list$down), "\n")
    ```

-  Which tissue shows more significantly accessible regions?
-  Is this in agreement with what you have seen in the Volcano plot?

## 4. Annotate Genomic Overlap with ChIPseeker

Next, you will annotate peaks to determine their genomic context (promoters, exons, introns, etc.):

**Annotation Parameters:**.  
- TSS region: ±1000 bp around transcription start sites.  
- Include detailed genomic annotations.  
- Separate analysis for up and down regulated peaks. 


```r
# Annotate all peaks
peakAnno = annotatePeak(gr, 
                        tssRegion=c(-1000, 1000), 
                        TxDb=TxDb, 
                        annoDb=AnnoDb, 
                        overlap = "TSS")


# Annotate upregulated peaks
peakAnno_up <- annotatePeak(gr_list$up, 
                           tssRegion = c(-1000, 1000), 
                           TxDb = TxDb, 
                           annoDb = AnnoDb, 
                           overlap = "TSS")

# Annotate downregulated peaks
peakAnno_down <- annotatePeak(gr_list$down, 
                             tssRegion = c(-1000, 1000), 
                             TxDb = TxDb, 
                             annoDb = AnnoDb, 
                             overlap = "TSS")

                             
# Save the objects
saveRDS(peakAnno, "results/06_DA_analysis/Annotated_peaks.rds")
saveRDS(peakAnno_up, "results/06_DA_analysis/Annotated_peaks_up.rds")
saveRDS(peakAnno_down, "results/06_DA_analysis/Annotated_peaks_down.rds")
```

**Task 3: Examine Annotation Results**

Explore the new objects created after annotation.   

- What is the frequency of promoters overlap overall?
- Do upregulated and downregulated peaks have different frequencies?  

 `peakAnno`, you can find different layers of information if you do `peakAnno@<tab>`:

```r
# Look at all peaks annotation summary
peakAnno
# Look at the upregulated peaks annotation summary
peakAnno_up
# Look at the downregulated peaks annotation summary
peakAnno_down
```

You can examine `peakAnno` object by doing `peakAnno@<tab>`, you will find several layers of information.  

- How many downregulated peaks overlap exons?  

```r
# Count how many downregulated peaks overlap exonic regions
sum(peakAnno_down@detailGenomicAnnotation$Exon)
```


**Task 4: Visualise Peak Annotations**


Use `plotAnnoPie` function from ChipSeeker package to visualize the genomic overlap distribution of upregulated peaks and downregulated peaks, separately:

??? Success "Solution"
    ```r
    # Plot pie charts for peak categories
    plotAnnoPie(peakAnno_up, main = "\n\nDA up") 
    plotAnnoPie(peakAnno_down, main = "\n\nDA down")
    ```

These plots show the percentage of peaks falling into different genomic categories (promoters, exons, introns, intergenic regions, etc.).

## 5. Functional Enrichment Analysis

Next, we will focus in peaks overlaping promoter regions (TSS) of genes, and we will perform functional enrichment analysis on those genes to understand which biological processes or metabolic pathways may be afected by the changes in chromatin accessibility.

**Task 6: Prepare gene list**

Start testing DA peaks that are upregulated (higher accessibility in Kidney).  

- Take DA upregulated peaks (peakAnno_up), select those overlapping TSS (Promoter), and get the list of gene SYMBOLs of the overlapping genes.  
- Get a list of universe gene SYMBOLs to compare the DA upregulated genes to.  

??? Success "Solution"
    ```r
    # Extract gene symbols for peaks in promoter regions
    genes_tss_up <- peakAnno_up@anno[peakAnno_up@detailGenomicAnnotation$Promoter, "SYMBOL"]

    # Create universe of genes with all promoter-overlapping peaks for background
    universe <- peakAnno@anno[peakAnno@detailGenomicAnnotation$Promoter,]$SYMBOL
    ```



**Task 7: Gene Ontology (GO) Enrichment Analysis**

Use `enrichGO` function from [`ClusterProfiler`](https://bioconductor.org/packages/devel/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf) to perform GO Enrichment Analysis of `genes_tss_up` set.   

- Provide unique universe. 
- Choose Biological Process (BP) ontology.  
- Apply multiple testing correction, and pvalueCutoff 0.05
- Pass AnnoDb variable created at the beginning of this tutorial section for OrgDb.  
- Specify keyType (SYMBOL).  

??? Success "Solution"
    ```r
    # GO enrichment for upregulated genes
    ego_up <- enrichGO(gene          = genes_tss_up$SYMBOL,
                    universe      = unique(universe),
                    OrgDb         = AnnoDb,
                    keyType       = "SYMBOL",
                    ont           = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

    # Check results
    ego_up
    ```


**Task 8**  

- Did you find significantly enriched GO terms?  

- Try to run the same analysis considering all genes as a universe, why do you think there is a difference?  

- Visualise the GO terms using a function from [ClusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) such as: `barplot()`

??? Success "Solution"
    ```r
    barplot(ego_up, showCategory = 20)
    ```

!!! Tip "Bonus"
    Repeat tasts 6-8 for DA downregulated peaks (higher accessiblility in Cerebrum)

    ??? Success "Solution"
        ```r
        # Create list of gene SYMBOLS associated to TSS overlapping peaks
        genes_tss_down <- peakAnno_down@anno[peakAnno_down@detailGenomicAnnotation$Promoter, "SYMBOL"]

        # Run enrichGO using gene SYMBOLS associated to all peaks overlapping promoters
        ego_down <- enrichGO(gene   = genes_tss_down$SYMBOL,
                    universe      = unique(universe),
                    OrgDb         = AnnoDb,
                    keyType       = "SYMBOL",
                    ont           = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

        # Plot results
        barplot(ego_down, showCategory = 20)
        ```




!!! tip "Next Steps"
    You can change the `ont` parameter to explore different ontologies:  
    - "MF": Molecular Function.  
    - "BP": Biological Process   
    - "CC": Cellular Component  


## 7. Additional Visualizations

There are several other plots that ChipSeeker allowes you to do. 
You can explore a couple of examples related to plotting the profile of peaks around certain regions


**Task 9: Peak Heatmap**

- Create a heatmap showing peak signal around TSS:  

```r
# Generate peak heatmap around gene bodies
peakHeatmap(peak = gr,
            TxDb = TxDb,
            upstream = 1000,
            downstream = 1000,
            by = "gene",
            type = "start_site",
            nbin = 800)
```

!!! tip "Bonus Task"

    Now let's plot the profile of peak counts around annotated TSS split by conditions

    ```r
    # Get TSS (promoter) annotation for mouse mm10
    promoter <- getPromoters(TxDb=TxDb, upstream=1000, downstream=1000)

    # Build a matrix of peaks per conditions, and combine them in a list
    tagMatrix_Cer <- getTagMatrix(gr_list$down , windows=promoter)
    tagMatrix_Kid <- getTagMatrix(gr_list$up , windows=promoter)
    tagMatrixList <- list(Cerebrum=agMatrix_Cer , Kidney=tagMatrix_Kid )

    # plot 
    plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")
    ```