# Clearwater
A web-deployed pipeline to analyze Legionella pneumophila samples.

# Workflow

```mermaid
gitGraph
       branch Clearwater_wf
       commit id: "Sequence fasta"

       branch MLST
       checkout MLST
       commit id: "elgato4mlst"
       commit id: "legsta"
       checkout Clearwater_wf
       merge MLST
       commit id: "Taxonomic classification" tag:"Mash"
       commit id: "Subspecies Classification" tag:"Kraken"
       branch phylogeny
       checkout phylogeny
       commit id: "MSA" tag:"ClustalW"
       commit id: "Phylogenetic tree" tag: "Iqtree"
       checkout Clearwater _wf
       merge phylogeny        
```

# Software Tools implemented/Dependencies
kraken, elgato, ksnp4, legsta, mash.

# Installation
### Clone repository
```
git clone https://github.com/BPHL-Molecular/clearwater.git
```

### Create your conda environment
