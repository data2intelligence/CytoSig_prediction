Example code of run CytoSig analysis and reproduce prediction results on bulk and single-cell cohorts  

Please CD into the src folder for the following steps.  

# Stage 1: download datasets  
run "./download.py"  
  
The output will be available in the folder data, named with bulk, single_cell, and output  

# Stage 2: convert single-cell data to CytoSig input format  
run "./convert_sc_data.py"  
**Note**: This step needs a CPU with large memory of 64G as the dataset EGAS00001004571 contains many single cells.  

# Stage 3: predict CytoSig signaling activity and generate figures  
run "./run.py"  
This step needs a CPU with large memory of 64G as the dataset EGAS00001004571 contains many single cells.

**Task 1**: bulk data from tumor and inflammatory disease studies    

**Task 2**: single-cell data from COVID19 studies    