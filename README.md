Example code of run CytoSig analysis and reproduce prediction results on bulk and single-cell cohorts  

Please CD into the src folder for the following steps.  

# Stage 1: download datasets  
run "./download.py"  
  
The starting data will be available in the folder data, named with bulk, single_cell, and output  

# Stage 2: convert single-cell data to CytoSig input format  
run "./convert_sc_data.py"  
**Note**: This step needs a CPU with large memory of 64G as the dataset EGAS00001004571 contains many single cells.  

We included two single-cell datasets that are neither in CytoSig input format or CellRanger format. GSE145926 is released as H5 files. EGAS00001004571 is released as Seurat object.  
This program will convert both datasets to python pickles of dense matrices as CytoSig input.  

# Stage 3: predict CytoSig signaling activity and generate figures  
run "./run.py"  
This step needs a CPU with large memory of 64G as the dataset EGAS00001004571 contains many single cells.  
The relevant figure numbers are labeled above each function.  

**Task 1**: bulk data from tumor and inflammatory disease studies  

**Task 2**: single-cell data from COVID19 studies  
