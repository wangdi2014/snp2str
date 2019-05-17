SNP2Str : Mapping PDB seq and UniProt seq 
=========================================

Authors: Lin An; Difei Wang, Ph.D.

Emails: difei.wang AT georgetown.edu or wang.difei AT yahoo.com or difei.wang AT nih.gov

Description
===========

The python code here was developed by Mr. Lin An and Difei Wang at Georgetown University.
This is the key part of our snp2str analysis package. If you use this code for your work, please cite the following paper. Thanks!

Wang, D. et al. SNP2Structure: A public and versatile resource for mapping and modeling msSNPs on human protein structures. Comput. and Struct. Biotech. J. 2015, 13, 514-519.	

Requirement
===========

  - Python and BioPython
  - Clustalw2

Installation
============
  - install miniconda or anaconda

  - conda create snp2str python=2.7
  - conda activate snp2str
  - conda install -c bioconda biopython
  - conda install -c bioconda clustalw

or 
  - conda create -f snp2str.yml

Running steps
=============
  - conda activate snp2str
  - cd mapping
  - ./run.sh
