#!/usr/bin/env python
import subprocess
import pandas as pd
import re
import os
import sys


def run_tblastn(dbpath, inputfile, outfile, evalue = "0.001", cpus = 1):
     out = subprocess.run(["tblastn", "-query", dbpath, "-subject", inputfile, "-out", outfile,
                           "-evalue", evalue, "-outfmt", "6", "-num_threads", str(cpus)])
     return(out)

def is_colocal(df, gene_class, int_sstart, trim_threshold):
    df_gene = df[df["Gene_class"] == gene_class]
    try:
        ret = abs(df_gene['sstart'][df_gene['evalue'].idxmin()] - int_sstart) < trim_threshold
    except:
        return(False)
    return(ret)

def gene_best_hit_location(df, gene):
    df_gene = df[df["Gene_class"] == gene]
    try:
        ret = df_gene['sstart'][df_gene['evalue'].idxmin()]
    except:
        return(0)
    return(ret)

def pici_direction(int_location, alpA_location):
     if alpA_location > int_location:
          return(1)
     else:
          return(-1)
     
def search(x):
    df = tblastn_results[tblastn_results['sseqid'] == x['sseqid']]
    direction = pici_direction(x['sstart'], gene_best_hit_location(df, "alpA"))
    result = {'Integrase': x['qseqid'], 'Integrase_pident': x['pident'],
              'int_location': x['sstart'],
              'alpA': is_colocal(df, "alpA", x['sstart'], trim_threshold),
              'alpA_location' : gene_best_hit_location(df, "alpA"), 
              'prirep': is_colocal(df, "pri-rep", x['sstart'], trim_threshold),
              'prirep_location': gene_best_hit_location(df, "pri-rep"),
              'sis': is_colocal(df, 'sis', x['sstart'], trim_threshold),
              'Contig': x['sseqid'], "Genome": x["Genome_name"],
              'direction': pici_direction(x['sstart'], gene_best_hit_location(df, "alpA")),
              "PICI_start": gene_best_hit_location(df, "pri-rep") - direction*attr_trim, "PICI_end": gene_best_hit_location(df, "pri-rep") + direction*attl_trim}
    return(result)

def classify(df):
     pici_class = 'Putative PICI'
     if not df['alpA'] and not df['prirep']:
          pici_class = "No PICI"
     if df["alpA"] and df["sis"] and df["prirep"]:
          pici_class = "SaPi"
     if df["alpA"] and df["prirep"] and pici_class != "SaPi":
          pici_class = "EcCi"
     result = {"Class": pici_class,
               'Integrase': df['Integrase'],
               'Direction': df['direction'],
               'Contig': df['Contig'],
               "PICI_start": df["PICI_start"],
               "PICI_end": df["PICI_end"]
     }
     return(result)


def flatten(t):
    return [item for sublist in t for item in sublist]

def preprocess_blast(path, genome = "blast"):
    tblastn_result = pd.read_csv(path, sep = "\t",  names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                                                'sstart', 'send', 'evalue', 'bitscore'])
    tblastn_result["Gene_class"] = [re.sub('[0-9]{1,2}', '', x) for x in tblastn_result['qseqid']]
    tblastn_result["int_match"]  = (tblastn_result["Gene_class"] == "int") & (tblastn_result['pident'] > int_threshold)
    tblastn_result["Genome_name"] = genome
    # intmatches = tblastn_result[tblastn_result['int_match'] == True]
    # intranges = [range(sorted(x)[0], sorted(x)[1]) for x in zip(intmatches['sstart'], intmatches['send'])]
    # overlaps = []
    # for i_index, i in enumerate(intranges):
    #     for j_index, j in enumerate(intranges):
    #         if i_index == j_index:
    #             continue
    #         if len(set(i).intersection(j)) > 0:
    #             overlaps.append((i_index, j_index))
    #TODO remove overlapping int hits
    return(tblastn_result)

def blast_dir(path):
     result_list = []
     for index, genome in enumerate(os.listdir(path)[0:3]):
          perc = index/len(os.listdir(path_path))*100
          if perc % 1 == 0:
               print( "\r" + str(perc) )
               sys.stdout.flush()
               r = run_tblastn("databases/BLAST_protein_db.faa", genome_path+genome, "tblastn.tsv", cpus = cpus)

          result_list.append(preprocess_blast("tblastn.tsv", genome))
    
     tblastn_results = pd.concat(result_list)
     return(tblastn_results)
     
int_threshold = 40
pri_rep_ident_threshold = 50

trim_threshold = 30000
pri_rep_trim = 25000
attr_trim = 2000
attl_trim = 25000
#genome_path = "data/Alfred/"
genome_path = ""
blast_path = "data/RefSeq/blast_small.tsv"
force = True
output_dir = "slet_virome_sensitivity"
cpus = 4
try:
    os.mkdir("output_dir")
except FileExistsError:
    if not force:
        exit(1)
    print("Directory already exists")        

if genome_path !=  "":
     tblastn_results = blast_dir(genome_path)
else:
     tblastn_results = preprocess_blast(blast_path)

pici_results =  tblastn_results[tblastn_results['int_match'] == True].apply(search, axis = 1, result_type = 'expand')

pici_classifications = pici_results.apply(classify, axis = 1, result_type = 'expand')
pici_classifications = pici_classifications.drop_duplicates(subset=['Contig', 'PICI_start'])
pici_classifications = pici_classifications[pici_classifications['Class'] != "No PICI"]
    

