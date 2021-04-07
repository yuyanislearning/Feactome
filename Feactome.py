import pandas as pd
import numpy as np
import math
import os
import json
import argparse
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser(description='get reactome analysis through API')
parser.add_argument('--input_file', type=str, help='input gene list', dest = 'file')
parser.add_argument('--total_file', type=str, help='total gene list', dest = 'total_file')
parser.add_argument('--output_file', type=str, help='output', dest = 'output')
parser.add_argument('--dir', type=str, help='data dir', dest = 'dir')
parser.add_argument('--fdr', type=float, help='fdr threshold', dest = 'fdr')


args = parser.parse_args()

# request from Reactome
os.system('cd ' + args.dir + ';'\
    ' curl -H "Content-Type: text/plain" --data-binary @'+args.file+' -X POST '\
        '--url https://reactome.org/AnalysisService/identifiers/projection/ > reactome_input_res.json')

# process json
with open(args.dir + '/reactome_input_res.json') as f:
    input_res = json.load(f)


# request from Reactome
os.system('cd ' + args.dir + ';'\
    ' curl -H "Content-Type: text/plain" --data-binary @'+args.total_file+' -X POST '\
        '--url https://reactome.org/AnalysisService/identifiers/projection/ > reactome_total_res.json')

# process json
with open(args.dir + '/reactome_total_res.json') as f:
    total_res = json.load(f)

# count lines
input_gene_num = 0
with open(args.dir + '/' + args.file) as f:
	for line in f:
		input_gene_num+=1

input_gene_num-=1

total_gene_num = 0
with open(args.dir + '/' + args.total_file) as f:
	for line in f:
		total_gene_num+=1

total_gene_num-=1

# index the total pathway for quick query
d_path = {}
i = 0
for pathway in total_res['pathways']:
	d_path[pathway['stId']] = i
	i+=1

p_l = []
pathid_l = []
pathname_l = []
E_l = []
Etotal_l = []


# for each pathway in input, find corresponding pathway in total
for pathway in input_res['pathways']:
	ID = pathway['stId']
	total_i = d_path[ID]
	k = pathway['entities']['found']
	n = total_res['pathways'][total_i]['entities']['found']

	p = 1 - hypergeom.cdf(k, total_gene_num, n, input_gene_num )
	p_l.append(p)
	pathid_l.append(ID)
	pathname_l.append(pathway['name'])
	E_l.append(k)
	Etotal_l.append(n)

fdr_l = multipletests(p_l, args.fdr, method = 'fdr_bh')[1]

df = pd.DataFrame({
	'pathway ID':pd.Categorical(pathid_l),
	'pathway name':pd.Categorical(pathname_l),
	'#Entities found':np.array(E_l),
	'#Entities total':np.array(Etotal_l),
	'FDR':fdr_l
	})

df = df.sort_values(by='FDR', ascending=True)
df = df[df.FDR<=args.fdr]

df.to_csv(args.output)
