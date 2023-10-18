# CITEViz

<img src="vignettes/citeviz-logo-light-cropped.png" alt="drawing" width="200"/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/maxsonBraunLab/CITEViz/actions/workflows/rcmdcheck.yaml/badge.svg)](https://github.com/maxsonBraunLab/CITEViz/actions/workflows/rcmdcheck.yaml)
<!-- badges: end -->

## Background

The rapid advancement of new genomic sequencing technology has enabled the development of multi-omic single-cell sequencing assays. These assays profile multiple modalities in the same cell and can often yield new insights not revealed with a single modality. For example, CITE-Seq (Cellular Indexing of Transcriptomes and Epitopes by Sequencing) simultaneously profiles the single-cell RNA transcriptome and the surface protein expression. The extra dimension of surface protein markers can be used to further identify cell clusters – an essential step for downstream analyses and interpretation. Additionally, multi-dimensional datasets like CITE-Seq require nuanced visualization methods to accurately assess the data. To facilitate cell cluster classification and visualization in CITE-Seq, we developed CITEViz. 

## About CITEViz

CITEViz is a single-cell visualization platform with a custom module that replicates the interactive flow-cytometry gating workflow. With CITEViz, users can investigate CITE-Seq specific quality control (QC) metrics, view multi-omic co-expression feature plots, and classify cell clusters by iteratively gating on the abundance of cell surface markers. CITEViz was developed to make multi-modal single-cell analysis accessible to a wide variety of biologists, with the aim to discover new insights into their data and to facilitate novel hypothesis generation. 

## Intended Audience
CITEViz is intended to be used by biologists familiar with flow cytometry and who can either 1) perform basic single-cell analysis or 2) collaborate in a team with a computational scientists. The input data format for CITEViz is a pre-processed Seurat object, which requires a basic level of coding skills in R and the ability to follow public Seurat vignettes. Since the essential feature of CITEViz is its iterative filtering process, it is not intended for any data preprocessing or normalization. In our experience, we found CITEViz to be a uniquely collaborative tool that leverages the strength of both bench and computational scientists to explore and analyze data together. 

## Citation

If you have used our package for your publications, please cite our paper: [CITE-Viz: Replicating the Interactive Flow Cytometry Workflow in CITE-Seq](https://www.biorxiv.org/content/10.1101/2022.05.15.491411v1)

## Usage

Please use the following link to start gating on Seurat-processed CITE-Seq data:

https://maxsonbraunlab.github.io/CITEViz/