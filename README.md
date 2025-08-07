# Sieler2025__ZF_StressHistory

## Project Overview

This repository contains the complete analysis pipeline and data for a zebrafish (Danio rerio) experiment investigating the effects of historical stress exposure on host-parasite interactions, microbiome composition, and gene expression patterns. The study examines how prior exposure to temperature stress and antibiotic treatment influences subsequent responses to parasite infection.

## Experimental Design

The experiment employed a factorial design with the following factors:
- **Historical Stressors**: Control (no prior stress), Temperature stress, Antibiotic treatment, or Combined (Temperature + Antibiotics)
- **Parasite Exposure**: Exposed or Unexposed to parasites
- **Time Points**: Baseline (0 days) and post-exposure (14 days)

### Treatment Groups
- **Control**: No prior stressor, no parasite exposure
- **Parasite**: No prior stressor, parasite exposed
- **Temperature**: Prior temperature stress, no parasite exposure  
- **Temperature_Parasite**: Prior temperature stress, parasite exposed
- **Antibiotics**: Prior antibiotic treatment, no parasite exposure
- **Antibiotics_Parasite**: Prior antibiotic treatment, parasite exposed
- **Antibiotics_Temperature**: Prior combined stress, no parasite exposure
- **Antibiotics_Temperature_Parasite**: Prior combined stress, parasite exposed

## Repository Structure

```
Sieler2025__ZF_StressHistory/
├── Code/                          # Analysis scripts and functions
│   ├── Analysis/                  # Main analysis workflows
│   │   ├── 16S/                   # Microbiome analysis scripts
│   │   ├── TaxProfiler/           # Taxonomic profiling
│   │   └── Transcriptomics/       # Gene expression analysis
│   ├── Dada2/                     # 16S sequencing processing
│   └── Functions/                 # Custom R functions
├── Data/                          # All project data
│   ├── Input/                     # Processed input data
│   ├── Metadata/                  # Sample metadata and experimental design
│   ├── Raw/                       # Raw sequencing data
│   └── R_objects/                 # R data objects
├── Results/                       # Analysis outputs and visualizations
│   ├── DEGxDAT/                   # Differential gene expression analysis
│   ├── DiffAbund/                 # Differential abundance analysis
│   ├── DiffExpGene/               # Gene expression results
│   ├── Infection_Mortality/       # Infection and mortality data
│   ├── Microbiome/                # Microbiome analysis results
│   ├── NeutralModel/              # Neutral community model analysis
│   └── DATxMortality/             # Disease-associated taxa correlations
├── Manuscript/                    # Manuscript drafts and supplementary materials
├── Biblio/                        # Bibliography and references
└── starting_doc.Rmd              # Getting started template
```

## Key Analyses

### 1. Microbiome Analysis (`Results/Microbiome/`)
- 16S rRNA gene sequencing analysis using DADA2 pipeline
- Taxonomic profiling and community composition analysis
- Differential abundance testing between treatment groups
- Neutral community model analysis

### 2. Transcriptomics Analysis (`Results/DEGxDAT/`)
- RNA-seq differential gene expression analysis
- Integration of gene expression with microbiome data
- Historical contingency effects on gene expression
- Recovery patterns post-parasite exposure

### 3. Host-Parasite Interactions (`Results/Infection_Mortality/`)
- Parasite load quantification
- Host mortality analysis
- Disease progression metrics

### 4. Disease-Associated Taxa (`Results/DATxMortality/`)
- Correlation analysis between specific microbial taxa and mortality
- Identification of potential disease-associated taxa

## Getting Started

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)
- Required R packages (see individual analysis scripts for specific dependencies)

### Setup
1. Clone this repository
<!-- 2. Open the project in RStudio using the `.Rproj` file -->
3. Run `starting_doc.Rmd` to set up the analysis environment
4. Install required packages as prompted

### Running Analyses
1. **Environment Setup**: Start with `starting_doc.Rmd` to load data and set paths
2. **Microbiome Analysis**: Navigate to `Code/Analysis/16S/` for microbiome processing
3. **Transcriptomics**: Use scripts in `Code/Analysis/Transcriptomics/`
4. **Results**: View generated HTML reports in the `Results/` directory

## Data Files

### Metadata (`Data/Metadata/metadata.csv`)
Contains complete experimental design information including:
- Sample IDs and tank assignments
- Treatment groups and historical stressors
- Time points and parasite exposure status
- Morphological measurements (weight, length)
- Sex and mortality data

### Key Variables
- `Treatment`: Current treatment group
- `History`: Historical stressor exposure (0=Control, 1=Single stressor, 2=Combined)
- `Time`: Time point (0=Baseline, 14=Post-exposure)
- `Pathogen`: Parasite exposure status (0=Unexposed, 1=Exposed)
- `Temperature`: Temperature stress (0=Control, 1=Stressed)
- `Antibiotics`: Antibiotic treatment (0=Control, 1=Treated)

## Output Files

### HTML Reports
- `Results/Microbiome/microbiome.html`: Complete microbiome analysis
- `Results/DEGxDAT/DEGxDAT_Integrated.html`: Integrated gene expression analysis
- `Results/SanityCheck.Rmd`: Data quality and experimental design validation

### R Objects
- Processed phyloseq objects for microbiome data
- DESeq2 results for transcriptomics
- Integrated analysis results

## Citation

If you use this code or data in your research, please cite:
```
Sieler, M. (2025). Historical stress exposure influences host-parasite interactions 
in zebrafish: An integrated microbiome and transcriptomics approach. 
[Manuscript in preparation]
```

## Contact

For questions about this repository or the analyses contained within, please contact:
- **Michael Sieler**: 

## License

This project is licensed under the terms specified in `LICENSE.md`.

---

**Note**: This repository contains large data files and analysis outputs. Some files may be stored using Git LFS (Large File Storage) to manage repository size effectively.
