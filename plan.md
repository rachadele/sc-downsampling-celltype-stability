
My supervisor wants to perform a computational experiment to determine the minimum sequencing depth required for accurate cell type calling in single-cell RNA sequencing data. He wants me to downsample Smart-seq data to various depths and evaluate the stability of cell type assignments. 

I took some notes on how to design the experiment:

- At what depth do cell type calls stabilize?
- is average 10x read depth enough for accurate cell type calling?

- Downsample to X reads per cell twice  
  - Two fake datasets  
  - All the same cells, produce two datasets  
  - At what depth does cell typing degrade?  
  - One cell type will be assigned two different types at some point depending on read depth  
  - E.g. at 1 read per cell, the same cell will be 2 different types in each dummy dataset  
  - What is the minimum X number of reads when this stops happening

I need help finding some test data first. Maybe we can use the human or mouse smart-seq reference data i'm using for another project:

[human](https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f)
download link: https://datasets.cellxgene.cziscience.com/dfba29ea-d368-44c2-bb35-69d3e8f730ce.h5ad
[mouse](https://cellxgene.cziscience.com/collections/e3aa612b-0d7d-4d3f-bbea-b8972a74dd4b)
download link: https://datasets.cellxgene.cziscience.com/6ce804da-13a1-4270-a726-cbeea399f7f1.h5ad

for now, start with human data, and downsample to ~1000 cells. use the cell_type field for cell type assignments.