# Clearwater
<img src="https://github.com/user-attachments/assets/6eecc007-cd6b-4ef3-ba2c-9bde5eef66ab" alt="Clear Water" width="15%" />

## Clearwater is a scalable, automated **Nextflow DSL2 pipeline** for genomic characterization, quality control, phylogenetic mapping, and machine learning-driven subtyping of ***Legionella pneumophila*** assemblies.

## 📖 What Clearwater Does

This pipeline provides a unified "sequence-to-insight" workflow for raw bacterial FASTA assemblies. It processes strains simultaneously through a hybrid architecture of linear analytical steps and asynchronous multi-sample cohort steps:

* **In Silico Serotyping & Subtyping:** Leverages `ELGATO` and `Legsta` for automated sequence typing (MLST).
* **Taxonomic Assignment:** Uses `Kraken2` for deep subspecies-level identification.
* **Assembly Quality Control:** Runs asynchronous `CheckM` evaluations to assess completeness and contamination.
* **Core & Pan-Genomic Phylogeny:** 
  * Identifies reference-free core SNPs and builds bootstrapped trees using `kSNP4` and `IQ-TREE`.
  * Generates structural annotations (`Prokka`) and executes gene-family clustering (`Roary`) to map the pangenome matrix.
* **Comparative Genomics:** Calculates Average Nucleotide Identity (`FastANI`) against reference libraries.
* **Predictive Analytics & Advanced Modeling:** 
  * Imputes missing loci elements (`miceforest`) and trains a **K-Nearest Neighbors (KNN)** classifier to resolve *L. pneumophila* subspecies.
  * Conducts **Principal Component Analysis (PCA)** post-classification to visualize high-dimensional feature distributions and variance.
  * Subjects downstream profiles to data-driven **Firth Logistic Regression** alongside a terminal **Lasso Multinomial Logistic Regression** model for feature pruning and classification.

---
## 📊 Pipeline Workflow Diagram
This chart outlines how independent single-sample data streams run alongside asynchronous, multi-sample cohort configurations:

```mermaid
graph TD
    %% Configuration and Styling
    classDef params fill:#f9f,stroke:#333,stroke-width:2px;
    classDef process fill:#bbf,stroke:#333,stroke-width:1px;
    classDef channel fill:#dfd,stroke:#333,stroke-width:1px;
    classDef collection fill:#ffd27f,stroke:#333,stroke-width:1px;

    %% Input Configurations / Parameters
    subgraph Input_Parameters [Pipeline Inputs]
        P1[params.input] -->|*.fasta| ch_samples(ch_samples channel)
        P2[params.ref_genomes_dir] --> ch_ref_dir(ch_ref_dir Value Channel)
        P3[params.krakenSub_db] --> kraken_db_file[krakenSub_db File]
    end

    %% Parallel Sample Processing Stream (Linear Chain)
    subgraph Sample_Stream [Per-Sample Parallel Processing Line]
        ch_samples --> elgato4mlst[elgato4mlst]
        elgato4mlst --> legsta[legsta]
        legsta -->|ch_legsta_out| kraken_subs[kraken_subs]
        kraken_db_file --> kraken_subs
        
        kraken_subs -->|ch_kraken_out| proc_mash[proc_mash]
        proc_mash --> pyProc1[pyProc1]
        pyProc1 --> prokka[prokka]
    end

    %% Asynchronous Batch / Collection Parallel Pathways
    subgraph Core_SNP_Phylogeny [Core Genome Assembly Branch]
        ch_samples -->|collect| ch_collected_fastas(Collected FASTAs)
        ch_collected_fastas --> ksnp4[ksnp4]
        ksnp4 -->|run_dir| iqtree_bootstrap[iqtree_bootstrap]
    end

    subgraph Assembly_QC [Asynchronous Quality Control]
        ch_collected_fastas --> checkm_qc[checkm_qc]
    end

    %% Gene Pangenome with Plotting Functionality
    subgraph Gene_Pangenome [Pangenome Extraction & Plotting Branch]
        prokka -->|ch_prokka_out| map_gff{Isolate GFF Files}
        map_gff -->|collect| ch_collected_gffs(Collected GFFs Array)
        ch_collected_gffs --> PANGENOME_ANALYSIS[PANGENOME_ANALYSIS]
        PANGENOME_ANALYSIS -->|run_dir| PLOT_PANGENOME[PLOT_PANGENOME]
    end

    %% Linear Chain Post-Prokka
    subgraph Individual_Annotation_and_ANI [ANI Verification Loop]
        prokka -->|ch_prokka_out| pyProc4[pyProc4]
        pyProc4 --> pyProc5[pyProc5]
        
        %% Combined Logic
        pyProc5 -->|ch_processed_samples| combine_checkm{combine}
        checkm_qc -->|results_dir| combine_checkm
        
        combine_checkm -->|ch_with_checkm| pyProc3[pyProc3]
        pyProc3 -->|ch_ani_inputs| runFastANI[runFastANI]
        ch_ref_dir --> runFastANI
    end

    %% Consolidation and Statistical Analytics Node
    subgraph Final_Reporting_Analytics [Reporting & Downstream Analytics]
        runFastANI -->|ch_ani_outputs| pyProc6[pyProc6]
        pyProc6 --> pyProc7[pyProc7]
        
        pyProc7 -->|ch_global_report| map_report{Map report.txt}
        map_report -->|collect| all_reports_ch(All Sample Reports Array)
        
        all_reports_ch --> compile_report[compile_report]
        compile_report -->|global_report| data_preprop[data_preprop]
        data_preprop -->|imputed_report| knn_classifier[knn_classifier]
        
        %% NEW STEPS ADDED HERE
        knn_classifier -->|knn_results| pca_analysis[pca_analysis]
        pca_analysis --> firth_regression[firth_regression]
        firth_regression --> SKLEARN_LASSO[SKLEARN_LASSO]
    end

    %% Assigning Classes for Visual Identification
    class P1,P2,P3 params;
    class elgato4mlst,legsta,kraken_subs,proc_mash,pyProc1,prokka,pyProc4,pyProc5,ksnp4,iqtree_bootstrap,checkm_qc,PANGENOME_ANALYSIS,PLOT_PANGENOME,pyProc3,runFastANI,pyProc6,pyProc7,compile_report,data_preprop,knn_classifier,pca_analysis,firth_regression,SKLEARN_LASSO process;
    class ch_samples,ch_ref_dir,ch_collected_fastas,ch_collected_gffs,all_reports_ch channel;
    class combine_checkm,map_gff,map_report collection;
```
---

## 🛠 Dependencies

The pipeline relies on a containerized or environment-managed runtime strategy to maximize reproducibility. You do not need to install the individual bioinformatics utilities manually.

### Core Frameworks
* **Nextflow** (Runtime tool; Java-based engine)

### Container Registries (Docker / Singularity)
The pipeline automatically fetches and launches isolated tool tasks using the following community layers:
* `staphb/prokka:latest` — Structural genome annotations
* `staphb/roary:latest` — High-performance pangenome calculation & graphing scripts
* `staphb/ksnp4:latest` — Core k-mer SNP phylogenetics and file layout parsing
* `staphb/kraken2:latest` — High-speed database string matching

### Conda Profiles (Local Fallbacks)
Statistical and machine learning scripts drop down to local environment recipes utilizing:
* `pandas` & `numpy` (Structured data parsing matrices)
* `scikit-learn` (KNN Classification engines, PCA dimensionality reduction, and Lasso Multinomial Logistic Regression modeling)
* `miceforest` (Iterative Chained Equations imputation)
* `matplotlib` & `seaborn` (Statistical matrix plotters & scatter/heatmap visualization engines)

---

## 💻 Resources Required

Requirements scale dynamically depending on your cohort batch sizing. Minimum allocations for the **10-sample testing baseline** include:

| Component | Minimum Specification | Recommended Specification | Notes |
| :--- | :--- | :--- | :--- |
| **CPU Threads** | 4 Cores | 8–16 Cores | `Roary (MAFFT)` and `kSNP4` scale linearly with core counts |
| **RAM** | 16 GB | 32 GB+ | Large `Kraken2` databases require memory headroom |
| **Storage** | 20 GB free space | 100 GB+ SSD | Accounting for database structures and intermediate caches |
| **OS** | Linux / macOS | Linux (Ubuntu/CentOS) | Essential for native Docker/Singularity daemon linking |

---

## 🚀 How to Run Clearwater(Usage)?

### 1. Clone the repository
```
git clone https://github.com/BPHL-Molecular/clearwater.git
```

### 2. Configure Parameters
Set up your file paths and parameters inside a params.yaml file or declare them inline during runtime execution.

### 3. Automatically create the pipeline conda environment to upload all conda packages
#### a) Execute this command on your terminal
```bash
conda env create -f CLEARWATERenvironment.yml
```
#### b) Activate the created environment
```bash
conda activate CLEARWATER
```
### 4. Standard Execution Command
Launch the master execution run by pointing to your local directory architecture or use the provided sbatch script( sbatch clearwater_run.sh):

```bash
nextflow run clearwater_wf.nf \
  --input "/path/to/your/fasta_folder" \
  --ref_genomes_dir "/path/to/reference_library" \
  --krakenSub_db "/path/to/kraken_database_file" \
  --output "results_output" \
  -profile apptainer
```

*To resume an interrupted or cached pipeline run without reprocessing completed steps, add the `-resume` flag.*

---

## 📂 Example of Output

When execution finishes, the designated output directory (`params.output/`) will be populated with structured computational folders:

```text
results_output/
├── [SampleID]/                          # Individual sample analytics folders
│   └── prokka_out/                      # GFF, FAA, FNA, and GBK file maps
├── pangenome_analysis/
│   └── roary_out/
│       ├── gene_presence_absence.csv    # Global pangenome matrix spreadsheet
│       ├── summary_statistics.txt       # Breakdown of Core, Shell, and Cloud distributions
│       └── core_gene_alignment.aln      # MAFFT core gene sequence alignment file
├── pangenome_plots/                     # Visualizations of the pangenome structure
│   ├── pangenome_frequency.png          # Gene frequency distribution plot
│   ├── pangenome_matrix.png             # Presence/absence matrix heatmap
│   └── pangenome_pie.png                # Core vs. Cloud distribution pie chart
├── pca_plots/                           # Dimensionality reduction outputs
│   └── subspecies_pca_comparison_panel  # Comparative PCA scatter layout across cohorts
├── knn_classification_results/          # Supervised machine learning performance data
│   ├── knn_performance_metrics.txt      # Model precision, recall, and accuracy scores
│   └── knn_predicted_subspecies         # Assigned taxonomic classifications per sample
├── lassoMregression_results/
│   ├── lasso_model_summary.txt          # Model feature selection stats
│   ├── lasso_coefficients_reportingTable.csv # Feature performance scores
│   └── lasso_coefficients_heatmap.png   # Multiclass predictive heatmap
├── firth_multinomial_formatted.txt      # Pre-formatted regression input data array
├── firth_multinomial_results.txt        # Firth-corrected logistic regression p-values and log-odds
├── ksnp_out/                            # Core k-mer SNP phylogenetics and validation files
└── bootstrap_tree/                      # IQ-TREE consensus tree and support value mappings
    ├── core_bootstrap.iqtree            # Detailed IQ-TREE execution log, models, and text-based tree
    └── core_bootstrap.treefile          # Maximum-likelihood consensus tree file (Newick format)
```
### 📈 Generated Visualizations & Performance

<details>
<summary><b>🔍 Click to view PCA Dimensionality Reduction Panel (pca_plots/)</b></summary>
<br>

To view the separation performance across *Legionella pneumophila* 7-gene multi-locus profiles, open this comparison panel:

<div align="center">
  <img src="https://github.com/user-attachments/assets/30067913-9d18-419b-80ed-a9cdc7a0240a" alt="subspecies_pca_comparison_panel" width="85%" />
  <p><i>Figure: Dimensionality reduction comparing full cohort ground truth (Left) against unseen kNN validation models (Right).</i></p>
</div>
</details>


<details>
<summary><b>📊 Click to view KNN Performance Metrics (knn_performance_metrics.txt)</b></summary>
<br>

Below is the evaluation summary and classification matrix generated by the machine learning module:

```text
======================================================================
         kNN SUBSPECIES CLASSIFIER PERFORMANCE SUMMARY                
======================================================================

Partitioning Strategy : 70% Training / 30% Testing
Overall Accuracy      : 66.67%

--- DETAILED EVALUATION METRICS (TEST SET) ---
                    precision    recall  f1-score   support

    subsp. fraseri       0.00      0.00      0.00         1
subsp. pneumophila       0.67      1.00      0.80         2

          accuracy                           0.67         3
         macro avg       0.33      0.50      0.40         3
      weighted avg       0.44      0.67      0.53         3

--- ALIGNED CONFUSION MATRIX ---
True / Pred        | subsp. fraseri     | subsp. pneumophila
------------------------------------------------------------
subsp. fraseri     | 0                  | 1                 
subsp. pneumophila | 0                  | 2                 
======================================================================
```
</details>

<details>
<summary><b>🧬 Click to view Pangenome Analysis Panel (pangenome_plots/)</b></summary>
<br>

This structured multi-panel visualization provides an overview of the global gene presence/absence patterns, frequency profiles, and pan-genome saturation distributions computed by the Roary analysis module:

<div align="center">
  <table border="0" cellspacing="0" cellpadding="5" style="border-collapse: collapse; border: none; width: 100%;">
    <tr style="border: none;">
      <td align="center" valign="top" style="border: none; width: 33.33%;">
        <b>Gene Frequency Distribution</b><br><br>
        <img src="https://github.com/user-attachments/assets/4d5076b2-53fe-4721-bb24-75290f9284b3" alt="pangenome_frequency" style="width: 100%; max-width: 350px; border-radius: 4px;" />
      </td>
      <td align="center" valign="top" style="border: none; width: 33.33%;">
        <b>Roary Alignment Matrix Map</b><br><br>
        <img src="https://github.com/user-attachments/assets/56554bd8-663f-414e-bab1-2daef4aef1d7" alt="pangenome_matrix" style="width: 100%; max-width: 350px; border-radius: 4px;" />
      </td>
      <td align="center" valign="top" style="border: none; width: 33.33%;">
        <b>Core vs. Cloud Distribution Pie</b><br><br>
        <img src="https://github.com/user-attachments/assets/341c1120-edb0-4445-a985-9398d3aa0ffa" alt="pangenome_pie" style="width: 100%; max-width: 350px; border-radius: 4px;" />
      </td>
    </tr>
  </table>
  <br>
  <p><i>Figure: Integrated pangenome analysis panel capturing core genome conservation metrics (1,280 core gene clusters), accessory shell trends (3,640 genes), and strain-specific cloud components (2,547 unique isolates) evaluated across the 10-strain cohort.</i></p>
</div>
</details>

<details>
<summary><b>🎯 Click to view LASSO Multinomial Regression Heatmap (lassoMregression_results/)</b></summary>
<br>

This heatmap maps the predictive weights (Log-Odds Coefficients) calculated by the machine learning module to identify significant diagnostic allele predictors for each target subspecies class:

<div align="center">
  <!-- Paste your exact unique asset link inside the src quotes below -->
  <img src="https://github.com/user-attachments/assets/e505bb46-ca2f-44f0-9828-3245027fe0c7" alt="LASSO Multinomial Regression Diagnostic Allele Predictors" width="75%" />
  <p><i>Figure: LASSO multiclass predictive heatmap tracking diagnostic genetic variants (e.g., Gene_neuA_neuAh: 2.0 strongly predicting subsp. pneumophila with a -2.8307 log-odds coefficient shift).</i></p>
</div>
</details>

<details>
<summary><b>🌳 Click to view Core Genome Phylogeny (bootstrap_tree/)</b></summary>
<br>

This Maximum-Likelihood consensus tree resolves the phylogenetic relationships across the target cohort based on core genome k-mer SNP alignments. High-confidence bootstrap support configurations (value of 100) are mapped onto the nodes:

<div align="center">
  <!-- Paste your new tree asset link inside the src quotes below -->
  <img src="https://github.com/user-attachments/assets/6edec3e2-150d-4bd6-b41e-9f810655e1fa" alt="Core Genome Phylogeny Tree" width="100%" style="background-color: white !important; padding: 15px; border-radius: 6px; border: 1px solid #e1e4e8;" alt="Core Genome Phylogeny Tree with Bootstraps" width="100%" style="background-color: white !important; padding: 15px; border-radius: 6px; border: 1px solid #e1e4e8;" />
  <br><br>
  
  <p><i>Figure: Midpoint-rooted phylogenetic tree computed via IQ-TREE. Branch text numbers indicate standard bootstrap support metrics calculated across 1,000 replicate pairings, demonstrating maximum confidence values across all major nodes.</i></p>
</div>
</details>
