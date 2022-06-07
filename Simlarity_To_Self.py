####################################### imports ########################################################################################

# <editor-fold desc="Imports">
import pandas as pd
import sys
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
from Bio import SeqIO
import json
import subprocess
import time
from multiprocessing import Pool
import regex as re
# </editor-fold>

####################################### functions ########################################################################################

# <editor-fold desc="Functions">
def dict_to_json_txt(d, file_path):
    json.dump(d, open(file_path, 'w'))

def json_txt_to_dict(file_path):
    d = json.load(open(file_path))
    return d

def send_needle_run(Human_peptidom_path, ref_peptide,n,path_to_save=None,save_json=False):
    """
    gets human peptidome and fasta peptides and return global alignment result using Needle Wununch algorithm , the results can be saved as Json file
    ref_peptide - path to Fasta file with the peptides we want to compare to the human peptidome
    Human_peptidom_path - path to Fasta file to human peptidome
    path_to_save path to save the dictionary with the result
    n - level of similarity that we want to define as a threshold
    """
    seq_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_peptide, "fasta")} # parse fasta file to dictionary
    dic = {}
    for pep in seq_dict.values():
        command= "needle -asequence asis:" +pep +" -bsequence {path} -gapopen 100 -gapextend 10 -endweight Y -endopen 100 -endextend 10 -error N -warning N -sprotein -stdout Y -filter Y| grep -E -B 6 -A 3 'Identity:\s+[{n}0-9]' | grep -E 'Identity:|2:|Score'".format(path=Human_peptidom_path,n=n)
        out,string_out=subprocess.getstatusoutput(command)
        if out==0: #if there is stdout
            if pep not in dic:
                dic[pep] = [string_out]
            else:
                dic[pep].append(pep)

    if save_json:
        dict_to_json_txt(dic,path_to_save)
    return dic


def send_needle_run_one_pep(pep,Human_peptidom_path="/home/perr/Downloads/similaritytoself/partial_human9mers.fasta"):

    """get the Human peptidome and one peptide each time and return dict of peptides withe similarity
    this function gets one peptide each iteration in order to make it easier to use with multiprocessing """
    dic={}
    command= "needle -asequence asis:" +pep +" -bsequence {} -gapopen 100 -gapextend 10 -endweight Y -endopen 100 -endextend 10 -error N -warning N -sprotein -stdout Y -filter Y| grep -E -B 6 -A 3 'Identity:\s+[6-9]' | grep -E 'Identity:|2:|Score:|Similarity'".format(Human_peptidom_path)
    out,string_out=subprocess.getstatusoutput(command)
    if out==0: #if there is stdout
        if pep not in dic:
            dic[pep] = [string_out]
        else:
            dic[pep].append(pep)

    return dic

ref_peptides="/run/user/1003/gvfs/afp-volume:host=HERTZ-LAB-NAS.local,user=elinorpe,volume=Elinor/June Analysis/semi_strict_25_steps_samples.fasta"
seq_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_peptides, "fasta")} # parse fasta file to dictionary


def main_func(seq_dict,path_to_save=None,save_file=False):
    """main function to initiate the multiprocessing similarity to seld"""
    pool = Pool()  # Create a multiprocessing Pool
    t1 = time.perf_counter()

    nested_dic=pool.map(send_needle_run_one_pep,seq_dict.values())

    t2 = time.perf_counter()
    print(f'Finished in {t2 - t1} seconds')

    main_dic={} # convert the nested dict to regular dict
    for dic in nested_dic:
        for k in dic:
            main_dic[k]=dic[k]
    if save_file:
        dict_to_json_txt(main_dic,path_to_save)
    return main_dic

def df_creator(main_dic):
    """get the similarity output and return the results as a df"""
    dic_data={}
    dic_data["match_peps"]=[]
    dic_data["superbinder"]=[]
    dic_data["similarity"]=[]
    dic_data["identity"]=[]
    dic_data["score"]=[]


    for k in main_dic:
          matched_peps=re.findall("[A-Z]{9}_",main_dic[k][0])
          matched_similarity=re.findall("Similarity:\s+(\d{1})",main_dic[k][0])
          matched_identity=re.findall("Identity:\s+(\d{1})",main_dic[k][0])
          matched_score=re.findall("Score:\s+\-*\d+\.\d+",main_dic[k][0])

          if len(matched_peps)>1:
              dic_data["superbinder"].extend([k]*len(matched_peps))
          else:
            dic_data["superbinder"].append(k)

          for i in matched_peps:
            dic_data["match_peps"].append(i.rstrip("_"))
          for j in matched_similarity:
                dic_data["similarity"].append(j)
          for m in  matched_identity:
                dic_data["identity"].append(m)
          for s in matched_score:
                dic_data["score"].append(s.lstrip("Score:"))
    df=pd.DataFrame.from_dict(dic_data)
    df = df.astype({"similarity": float,"identity":float,"score":float})
    return df
# </editor-fold>






