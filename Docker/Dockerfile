FROM rocker/rstudio:4.5

RUN apt-get update && \
    apt-get install -y \
        libxt6 \
        libgl1 \
        libcurl4-openssl-dev \
        libglpk-dev \
        libxml2-dev \
        libproj-dev \
        libudunits2-dev \
        libgdal-dev \
        libcairo2-dev \
        libxt-dev \
        libz-dev \
        libbz2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add Miniconda to PATH for subsequent commands
ENV PATH="/opt/miniconda3/bin:${PATH}"
#ENV CONDA_ENV_NAME=atac_env

# Install miniconda and macs3 and the preredoquisites (cython and numpy)
RUN MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" && \
    wget $MINICONDA_URL -O miniconda.sh  && \
   # mkdir miniconda3 && \
   # mkdir -p /root/.conda && \
    bash miniconda.sh -b -p /opt/miniconda3 && \
    rm -f miniconda.sh && \
    conda clean -afy && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda create -y -n atac_env && \
    conda install -y -n atac_env -c bioconda -c conda-forge ataqv picard multiqc deeptools cython numpy samtools bedtools subread && \
    conda run -n atac_env pip install macs3==3.0.3 && \
    conda clean -afy

# Ensure conda is initialized for bash in the rstudio user's environment
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/rstudio/.bashrc 
#RUN conda init bash 

# Install R packages
RUN Rscript -e "install.packages('BiocManager'); BiocManager::install(c('GenomicAlignments','ggplot2','rtracklayer','dplyr','EnhancedVolcano','ChIPseeker', 'pheatmap','RColorBrewer', 'clusterProfiler','enrichplot'))"
RUN Rscript -e "BiocManager::install(c('DESeq2','apeglm','TxDb.Mmusculus.UCSC.mm10.knownGene','org.Mm.eg.db','rGREAT','GenomicRanges'))"
