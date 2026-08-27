# Spatial and TCR Profiling of GPH Plus Immune Checkpoint Therapy in Pancreatic Ductal Adenocarcinoma

## Overview

Pancreatic ductal adenocarcinoma (PDAC) is one of the most aggressive malignancies and is characterized by a highly immunosuppressive tumor microenvironment (TME). Our previous studies demonstrated that the combination of **gemcitabine, hydroxychloroquine, and paricalcitol (GPH)** can modulate the immune microenvironment in PDAC.

In this study, we investigated the therapeutic and immunological effects of combining **GPH with immune checkpoint therapy (IT; anti-PD-1 and anti-CTLA-4)**. We integrated **SpaTial Enhanced Resolution Omics-sequencing (Stereo-seq)** with **single-cell T-cell receptor (scTCR) sequencing** to characterize the spatial, cellular, and clonal remodeling of the PDAC tumor microenvironment following treatment.

## Study Design

Orthotopic PDAC models were evaluated across four treatment conditions:

* **Sham** – control treatment
* **IT** – anti-PD-1 + anti-CTLA-4 immune checkpoint therapy
* **GPH** – gemcitabine + hydroxychloroquine + paricalcitol
* **GPH+IT** – combined GPH and immune checkpoint therapy

Stereo-seq was used to characterize treatment-associated spatial changes within the tumor microenvironment, while single-cell TCR sequencing was used to investigate T-cell clonality, expansion, phenotype, and treatment-associated immune responses.

## Key Findings

The **GPH+IT combination produced the strongest antitumor response**, with greater tumor growth inhibition than the individual treatment groups. Survival studies further demonstrated the superiority of GPH+IT, which was the only treatment condition associated with complete tumor remission.

T cells isolated from treated mice retained antitumor activity and controlled orthotopically implanted PDAC tumors in both immunocompetent and RagKO mice.

Single-cell TCR analysis revealed preferential enrichment of T-cell clones, particularly **effector CD4+ T-cell clones**, within tumor–immune interface regions. These clones underwent substantial expansion following T-cell infusion. Rather than displaying predominantly terminal exhaustion, expanded clones showed transcriptional features associated with **effector activation**, together with increased **memory CD8+ T-cell** and **proliferative T-cell** populations.

Stereo-seq analysis demonstrated extensive spatial remodeling of the PDAC microenvironment following GPH+IT treatment. Major observations included:

* M1 macrophages increased from undetectable levels in sham-treated tumors to approximately **7.5% of all profiled cells**, representing approximately **83% of macrophages** following GPH+IT treatment.
* **CD8+ and cytotoxic T cells localized closer to tumor regions**, consistent with increased immune accessibility.
* Tumor-associated **hypoxia scores decreased** following combination treatment.
* Spatial ecotype analysis showed that GPH+IT reorganized tumors into more **spatially open and immune-accessible configurations** that were not observed in the other treatment conditions.

## Biological Significance

Together, these findings demonstrate that **GPH+IT induces coordinated clonal, cellular, and spatial remodeling of the PDAC tumor microenvironment**. The combination treatment promotes immune infiltration, shifts macrophages toward an M1-like state, reduces tumor-associated hypoxia, and supports the expansion and activation of antitumor T-cell clones.

Integration of Stereo-seq and single-cell TCR profiling provides a framework for connecting **T-cell clonal dynamics with their spatial organization and functional state within the tumor microenvironment**.

These findings provide mechanistic insight into the therapeutic activity of GPH combined with immune checkpoint blockade and may help identify spatial and immune biomarkers associated with treatment response. The results also provide a rationale for further investigation of this therapeutic strategy in combination with **TCR-based immunotherapies and adoptive T-cell transfer approaches**.

## Repository Contents

This repository contains the computational workflows and analysis code used for:

* Stereo-seq preprocessing and quality control
* Spatial cell-type characterization
* Tumor–immune spatial interaction analysis
* Spatial ecotype analysis
* Tumor hypoxia analysis
* Single-cell TCR repertoire analysis
* T-cell clonotype expansion analysis
* T-cell phenotype and functional-state characterization
* Integration of spatial transcriptomic and TCR-associated findings
* Generation of manuscript figures and supporting analyses

## Citation

Citation information will be updated upon publication of the associated manuscript.
