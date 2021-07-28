import pandas as pd
import numpy as np
import math
import os
import json
import argparse
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser(description='get reactome analysis through API')
parser.add_argument('--input_files', type=str, help='input gene list', dest = 'file')
parser.add_argument('--output_dir', type=str, help='output', dest = 'outdir')
parser.add_argument('--dir', type=str, help='data dir', dest = 'dir')
parser.add_argument('--fdr', type=float, help='fdr threshold', dest = 'fdr')

# example of json 
#>>> dat['pathways'][0]
#{'stId': 'R-HSA-381426', 'dbId': 381426, 'name': 'Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)', 'species': {'dbId': 48887, 'taxId': '9606', 'name': 'Homo sapiens'}, 'llp': True, 'entities': {'resource': 'TOTAL', 'total': 127, 'found': 6, 'ratio': 0.008733324164489065, 'pValue': 2.573869459165934e-05, 'fdr': 0.00622876409118156, 'exp': []}, 'reactions': {'resource': 'TOTAL', 'total': 14, 'found': 2, 'ratio': 0.0010358094110683633}, 'inDisease': False}

args = parser.parse_args()

ids = []
names = []
FDRs = []

# request from Reactome
with open(args.dir+'/'+args.file) as f:
    for file in f:
        os.system('cd ' + args.dir + ';'\
            ' curl -H "Content-Type: text/plain" --data-binary @\''+file+'\' -X POST '\
                '--url https://reactome.org/AnalysisService/identifiers/projection/ > temp.json')

        with open(args.dir + '/temp.json') as f:
            input_res = json.load(f)


        # for each pathway in input, find corresponding pathway in total
        for pathway in input_res['pathways']:
            ids.append(pathway['stId'])
            names.append(pathway['name'])
            FDRs.append(pathway['entities']['fdr'])


        df = pd.DataFrame({
            'pathway ID':pd.Categorical(ids),
            'pathway name':pd.Categorical(names),
            'FDR':FDRs
            })



        df = df.sort_values(by='FDR', ascending=True)
        df = df[df.FDR<=args.fdr]

        df.to_csv(os.path.join(args.dir,args.outdir,file[0:len(file)-4]+'.csv'))

