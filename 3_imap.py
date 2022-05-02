#### IMPORT LIBRARY ####
import scanpy as sc
import  imap
import pandas as pd
import numpy as np

adata = sc.read_loom('sdata220115.loom',sparse=False)
adata.obs[['ClusterName']].to_csv('cellid.csv',index= True ,header=False)  
adata.obs[['ClusterName']]

# add batch
adata.obs['batch']='4'
adata.obs['batch'][140:189]='1'
adata.obs['batch'][189:351]='2'
adata.obs['batch'][351:542]='3'

adata = imap.stage1.data_preprocess(adata, 'batch',n_top_genes=20000)  
# Preprocess the data.(count data -> log-format data, high-var genes selected)
# Output the basic information of the preprocessed data.


### Stage I
EC, ec_data = imap.stage1.iMAP_fast(adata, key="batch", n_epochs=150) 
### Stage II
output_results = imap.stage2.integrate_data(adata, ec_data, n_epochs=150)
output11=pd.DataFrame(output_results)
output11.columns=adata.var_names  
output11.index=adata.obs_names
output11.to_csv('output11.csv',index= True ,header=True) 
