#!/usr/bin/env python


#import libraries
import pandas as pd
import numpy as np
import io
import os
import sys 
import warnings

warnings.filterwarnings("ignore")

#now we read the reults of the logistic regression
results_logistic_regression_path=sys.argv[1]
results_logistic_regression=pd.read_table(results_logistic_regression_path, sep='\s+')
##now we will also read the freqx file that contains the allele frequencies
freqx_path=sys.argv[2]
freq_report=pd.read_table(freqx_path, sep='\t')
#we select only the columns that we need
freq_report=freq_report[['SNP','C(HOM A1)', 'C(HET)', 'C(HOM A2)',
       'C(HAP A1)', 'C(HAP A2)', 'C(MISSING)']]
#now we will merge the results with the variant info and the allele frequencies
variantInfo_results_and_freqx_merged=pd.merge(results_logistic_regression, freq_report, on='SNP', how='inner')
#no we want to create a new columns that contains the hyperlink to dbSNP
variantInfo_results_and_freqx_merged['dbSNP_link']='https://www.ncbi.nlm.nih.gov/snp/'+variantInfo_results_and_freqx_merged['SNP']
#now we store this in a csv
path_to_store=sys.argv[3]
variantInfo_results_and_freqx_merged.to_csv(path_to_store, index=False)

#now we will select only those values that are significant
significant_variants=variantInfo_results_and_freqx_merged[variantInfo_results_and_freqx_merged['P']<0.0000005]
#we will store this in a csv
significant_variants_path=sys.argv[4]
significant_variants.to_csv(significant_variants_path, index=False)

#we will also get the top ten most significant variants
topten=variantInfo_results_and_freqx_merged.nsmallest(10, 'P')
#we will store this in a csv
topten_path=sys.argv[5]
topten.to_csv(topten_path, index=False)
