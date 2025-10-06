# mechanism_to_model

## Overview

This repository contains bioinformatics analyses investigating the relationship between epigenetic mechanisms and physiological states in coral species. The project integrates multi-omics data including RNA-seq, whole genome bisulfite sequencing (WGBS), and physiological measurements to understand how epigenetic modifications correlate with coral stress responses and energetic states.

The primary focus is on two coral species:
- **Pocillopora spp.** - studied in the context of chronic low-level nutrient enrichment
- **Acropora pulchra** - studied in the context of simulated marine heatwave stress

## Key Features

- **Physiological Data Analysis**: Principal component analysis (PCA) of physiological variables to create proxies for energetic state
- **Epigenetic Annotations**: BLAST-based annotation of epigenetic machinery proteins using UniProt/Swiss-Prot databases
- **DNA Methylation Analysis**: Identification of differentially methylated loci (DML) using methylKit
- **Gene Expression Analysis**: Differential gene expression analysis and correlation with physiological states
- **GO Term Enrichment**: Gene ontology analysis for biological process interpretation
- **Integrative Analysis**: Correlation between epigenetic protein expression and physiological parameters

## Repository Structure

```
mechanism_to_model/
├── data/                           # Input data files
│   ├── 01-Physiology/             # Physiological measurements
│   ├── 03-HistonePTMs(simualtedData)/ # Histone post-translational modification data
│   ├── RNAseq/                    # RNA-seq data files
│   ├── *.csv                      # Various data matrices and metadata
│   └── *.faa                      # Protein sequence files
├── scripts/                        # Analysis scripts (R Markdown)
│   ├── 01-Physiology.Rmd          # Physiological PCA analysis
│   ├── 02-Pver-epimods-blast.Rmd  # Pocillopora annotation
│   ├── 03-EpigeneticProteins-Phys-Correlation.Rmd # Correlation analysis
│   ├── 04-carbon-transport.Rmd    # Carbon transport pathway analysis
│   ├── 05-methylKit.Rmd           # DNA methylation analysis
│   ├── 06-moreGO-annotations.Rmd  # Additional GO annotations
│   ├── 07-deg.Rmd                 # Differential expression analysis
│   └── 08-Apul-epimods-blast.Rmd  # Acropora annotation
├── output/                         # Analysis results and outputs
│   ├── 02-Pver-blast/             # BLAST results for Pocillopora
│   ├── 04-carbon-transport/       # Carbon transport analysis results
│   ├── 05-methylKit/              # Methylation analysis results
│   ├── 06-moreGO/                 # GO enrichment results
│   ├── 07-degs/                   # Differential expression results
│   └── 08-Apul-epimods-blast/     # BLAST results for Acropora
└── README.md                       # This file

```

## Data Description

### Input Data Sources

1. **Physiological Data**
   - Pocillopora: From [Becker et al. 2021](https://link.springer.com/article/10.1007/s00338-021-02138-2)
     - Repository: [Chronic_low_nutrient_enrichment](https://github.com/daniellembecker/Chronic_low_nutrient_enrichment_benefits_coral_thermal_performance_fore_reef_habitat)
   - Acropora pulchra: From Becker's 2022 heatwave experiment
     - Repository: [A.pul_Heatwave](https://github.com/daniellembecker/A.pul_Heatwave)

2. **Molecular Data**
   - RNA-seq gene count matrices
   - WGBS methylation data (Bismark output)
   - Protein sequences from genome annotations

3. **Reference Data**
   - UniProt/Swiss-Prot annotated protein database
   - Gene ontology (GO) term databases
   - Species-specific genome annotations

## Analysis Workflow

### 1. Physiological State Assessment (`01-Physiology.Rmd`)
- Calculate PC1 from physiological variables (biomass, chlorophyll, symbiont density, Pmax)
- Create proxy metrics for coral energetic state

### 2. Epigenetic Machinery Annotation (`02-Pver-epimods-blast.Rmd`, `08-Apul-epimods-blast.Rmd`)
- BLAST protein sequences against Swiss-Prot database
- Identify epigenetic machinery proteins (writers, readers, erasers)
- Extract GO term annotations

### 3. Correlation Analysis (`03-EpigeneticProteins-Phys-Correlation.Rmd`)
- Correlate epigenetic protein expression with physiological PC1
- Identify relationships between epigenetic regulation and energetic state

### 4. Pathway Analysis (`04-carbon-transport.Rmd`)
- Analyze specific biological pathways (e.g., carbon transport)
- Generate heatmaps of gene expression patterns

### 5. DNA Methylation Analysis (`05-methylKit.Rmd`)
- Process Bismark WGBS output files
- Identify differentially methylated loci (DML)
- Visualize methylation patterns and coverage
- Based on [Yaamini Venkataraman's methylKit workflow](https://osf.io/u46xj)

### 6. Gene Ontology Enrichment (`06-moreGO-annotations.Rmd`)
- Perform GO term enrichment analysis
- Interpret biological processes affected by treatments

### 7. Differential Gene Expression (`07-deg.Rmd`)
- Merge DEG data with GO annotations
- Generate heatmaps of differentially expressed genes
- Analyze treatment effects at the transcriptional level

## Dependencies

### Software Requirements
- R (≥ 4.0)
- NCBI BLAST+ (2.15.0+)
- Bismark (for WGBS data processing)

### R Packages

**Data manipulation and visualization:**
- tidyverse
- readr
- readxl
- reshape2
- knitr
- kableExtra
- DT

**Bioinformatics:**
- BiocManager
- Biostrings
- GenomicRanges
- methylKit
- genefilter

**Visualization:**
- ComplexHeatmap
- vegan
- gplots
- dichromat

**Statistical analysis:**
- stats4
- tm

### Installing R Dependencies

```r
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("methylKit", "Biostrings", "GenomicRanges"))

# Install CRAN packages
install.packages(c("tidyverse", "readxl", "vegan", "gplots", 
                   "dichromat", "ComplexHeatmap", "genefilter"))
```

## Getting Started

### Prerequisites
1. Clone this repository
2. Ensure you have R and required packages installed
3. Install NCBI BLAST+ if running annotation scripts
4. Download reference databases (Swiss-Prot, GO terms) if needed

### Running Analyses

The scripts are numbered to reflect the recommended execution order:

1. **Start with physiological data processing:**
   ```r
   # Open in RStudio
   rmarkdown::render("scripts/01-Physiology.Rmd")
   ```

2. **Run species-specific annotations:**
   ```r
   # For Pocillopora
   rmarkdown::render("scripts/02-Pver-epimods-blast.Rmd")
   
   # For Acropora
   rmarkdown::render("scripts/08-Apul-epimods-blast.Rmd")
   ```

3. **Continue with downstream analyses** in numerical order (03-07)

### Output Files

Each script generates outputs in the corresponding `output/` subdirectory:
- Annotated gene lists (`.tsv`, `.csv`)
- Figures (`.jpeg`, `.png`)
- Statistical results
- R data objects (`.RData`)

## Key Data Files

- `data/metadata.WGBS.csv` - Sample metadata for WGBS samples
- `data/DEGSeq2.sig.results.host.csv` - Significant differential expression results
- `data/annot_GO.terms.host.csv` - GO term annotations for host genes
- `data/merged_DEG_full_annot.csv` - Merged differential expression with full annotations
- `data/Pver_proteins_names_v1.0.faa` - Pocillopora protein sequences
- `data/Machinery.fasta` - Epigenetic machinery protein sequences

## Usage Tips

1. **File Paths**: Scripts use relative paths from the `scripts/` directory. Run scripts from within the `scripts/` folder or adjust paths accordingly.

2. **Evaluation Control**: Many code chunks in Rmd files have `eval = FALSE`. Set `eval = TRUE` to run specific analyses.

3. **Computational Resources**: BLAST and methylKit analyses can be computationally intensive. Adjust thread counts in BLAST commands (`-num_threads`) based on your system.

4. **Large Data Files**: WGBS data files are gitignored (see `.gitignore`). Contact the repository maintainers for access to raw sequencing data.

## Contributing

Contributions to improve the analyses or documentation are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/YourFeature`)
3. Commit your changes with clear messages
4. Push to your branch (`git push origin feature/YourFeature`)
5. Open a Pull Request

### Code Style
- Follow tidyverse style guidelines for R code
- Document code chunks with clear comments
- Use descriptive variable names
- Include figure captions and axis labels

## Citation

If you use this repository or build upon these analyses, please cite:

- Becker, D. M., et al. (2021). Chronic low level nutrient enrichment benefits coral thermal performance in fore reef habitat. *Coral Reefs*, 40(5), 1637-1655. [https://doi.org/10.1007/s00338-021-02138-2](https://doi.org/10.1007/s00338-021-02138-2)

Additional analyses in this repository are part of ongoing research. Please contact the repository maintainers for citation information.

## Related Resources

- [methylKit documentation](https://bioconductor.org/packages/release/bioc/html/methylKit.html)
- [Yaamini Venkataraman's methylKit workflow](https://osf.io/u46xj)
- [UniProt REST API documentation](https://www.uniprot.org/help/api)
- [NCBI BLAST+ documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## Contact

For questions or collaboration inquiries, please open an issue in this repository.

## License

Please refer to the repository license file for usage terms and conditions.

---

**Last Updated**: October 2024