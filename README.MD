# Slow vs Fast Codon Usage in Inter-Domain Connecting Regions  
### A Codon-Level Analysis of Multi-Domain Histidine Kinases in *Bacillus subtilis*

## 📌 Project Overview

This repository contains a complete, reproducible bioinformatics pipeline to investigate whether **slow codons are enriched in inter-domain connecting regions** of multi-domain proteins.

The analysis focuses on **20 sensor histidine kinases** from *Bacillus subtilis* subsp. *subtilis* strain 168, a classical model organism for bacterial signal transduction.

The central biological question is:

> **Are slow-translating codons preferentially located at inter-domain connecting regions, potentially facilitating co-translational folding?**

This hypothesis is tested using **codon-level annotation**, **domain-aware region mapping**, and **rigorous statistical validation**, including exact tests and permutation-based null models.

---

## 🧬 Biological Background

Synonymous codons are translated at different rates due to differences in tRNA abundance and codon usage bias. In some systems, clusters of slow codons have been proposed to induce ribosomal pausing, allowing nascent protein domains to fold before translation proceeds further.

Multi-domain proteins—such as histidine kinases—provide a natural framework to test this idea, as they contain clearly separated functional domains connected by linker regions.

---

## 🧪 Dataset

- **Organism:** *Bacillus subtilis* subsp. *subtilis* strain 168  
- **Protein family:** Sensor histidine kinases (two-component systems)  
- **Number of proteins:** 20  
- **Data sources:**
  - UniProtKB (protein sequences and domain annotations)
  - NCBI RefSeq / GenBank (coding sequences and genomic coordinates)

All coding sequences were retrieved using exact genomic coordinates (`coded_by` qualifiers) and **validated by exact translation match to UniProt protein sequences**.

---

## ⚙️ Analysis Pipeline

The pipeline consists of the following steps:

1. **CDS Retrieval and Validation**
   - UniProt → RefSeq cross-referencing
   - NCBI EFetch (`coded_by`) parsing
   - Strand-aware CDS extraction
   - Translation using bacterial genetic code (table 11)
   - Validation against UniProt protein sequences

2. **Domain Feature Extraction**
   - Domain-like annotations retrieved via UniProt REST API
   - Feature types include: Domain, Region, Transmembrane, Topological domain, Repeat, Motif, Coiled coil

3. **Definition of Connecting Regions**
   - Regions between adjacent domains
   - ±5 amino acid margin applied to account for boundary uncertainty
   - Mapping from amino-acid space to codon indices

4. **Codon Classification**
   - Codon adaptiveness weights (*w*) computed using a CAI-like framework
   - **Slow codons:** bottom 25% of *w* values
   - **Fast codons:** top 25% of *w* values

5. **Statistical Analysis**
   - 2×2 contingency tables (connecting vs non-connecting × slow vs not slow)
   - Odds ratio (effect size)
   - Fisher’s exact test (two-sided and one-sided)
   - 95% confidence interval for OR
   - **Within-gene permutation test (B = 5000)** to validate results under a null model

---

## 📊 Key Results

- Total codons analyzed: **9,735**
- Slow codon proportion:
  - Connecting regions: **15.99%**
  - Non-connecting regions: **16.06%**
- **Odds Ratio:** 0.995  
- **95% CI:** [0.869, 1.140]  
- **Fisher’s exact test:** p = 0.972 (two-sided)
- **Permutation test:** observed OR lies near the center of the null distribution

📌 **Conclusion:**  
Slow codons are **not enriched** in inter-domain connecting regions of *B. subtilis* histidine kinases.

---

## 📁 Repository Structure

```text
.
├── cds_from_ncbi.fna                      # Retrieved CDS sequences (FASTA)
├── cds_fetch_log.tsv                      # CDS retrieval and validation log
├── uniprot_domain_features.tsv            # Domain annotations from UniProt
├── per_codon_table.tsv                    # Codon-level annotations
│
├── download_uniprot_features.py           # UniProt feature retrieval
├── 02_get_cds_from_ncbi.py                # CDS retrieval from NCBI
├── 05_slow_fast_connecting_stats.py       # Core codon + region analysis
├── stats_fisher_ci.py                     # Fisher exact test + CI
├── permutation_test_or.py                 # Permutation-based validation
├── make_figures.py                        # Figure generation scripts
│
├── slow_fast_codon_connecting_report_final.docx
├── slow_fast_codon_connecting_report_SUPPLEMENTARY.docx
│
└── README.md
```
## 📈 Figures

The repository includes scripts to generate:

- Figure 1: Proportion of slow codons in connecting vs non-connecting regions

- Figure 2: Distribution of codon adaptiveness weights (w)

- Figure 3: Permutation-based null distribution of odds ratios

- Figures are embedded in the final Word report.

## 🔬 Reproducibility

- All data retrieval is performed programmatically using public APIs (UniProtKB, NCBI).

- All intermediate files are saved and versionable.

- The analysis can be reproduced end-to-end using the provided scripts.

## Requirements

- Python ≥ 3.9

- Biopython

- pandas

- numpy

- scipy

- statsmodels

- matplotlib

Install dependencies:

````bash
pip install biopython pandas numpy scipy statsmodels matplotlib
````

## 🎓 Scientific Significance

This project provides a well-controlled negative result, demonstrating that slow-codon enrichment at inter-domain boundaries is not a universal feature of bacterial multi-domain proteins.

The analysis highlights the importance of:

- effect-size interpretation,

- appropriate null models,

- and avoiding overgeneralization of co-translational folding hypotheses.

## 📚 References

Sharp, P. M., & Li, W.-H. (1987). The codon adaptation index. Nucleic Acids Research.

Jansen, R., Bussemaker, H. J., & Gerstein, M. (2003). Revisiting the CAI. Nucleic Acids Research.

The UniProt Consortium. UniProt: the universal protein knowledgebase.