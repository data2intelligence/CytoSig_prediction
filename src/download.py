#!/usr/bin/env python

import os, pathlib, tarfile
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

for url in ['https://hpc.nih.gov/~Jiang_Lab/CytoSig/CytoSig_prediction/data.tar.gz']:
    f = os.path.basename(url.rstrip('/'))
    
    out = os.path.join(base_path, f)
    urllib.request.urlretrieve(url, out)
    
    if url.find('.tar.gz') > 0:
        my_tar = tarfile.open(out)
        my_tar.extractall(base_path)
        my_tar.close()
        
        os.remove(out)
