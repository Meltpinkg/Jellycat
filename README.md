# Jellycat 
---	
### Installation 
	$ git clone https://github.com/Meltpinkg/Jellycat.git && cd Jellycat/

---	
### Introduction 

Jellycat is an ultra-fast and sensitive integration toolset especially for population-based SV calling and large callsets merging.
---
### Dependence 
	
	1. python3
	2. pysam
	3. numpy
---
### Usage 
	python cuteSV_merge.py <input.txt> <output.vcf> <work_dir>


| Parameter | Description | Default |
| :------------ |:---------------|-------------:|
|--threads|Number of threads to use.| 16 |
|--annotation| Annotation file to add.|NULL|
|--massive| Choice for large sample. |False|
|--support| Minimum support sample number to report. |1|
|-max_dist| Maximum distance for merging two SVs.|1000|

---
### Contact 
For advising, bug reporting and requiring help, please post on [Github Issue](https://github.com/Meltpinkg/Jellycat/issues) or contact sakura_qian@163.com.