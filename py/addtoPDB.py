import sys
import inspect
from glob import glob
import numpy
import os
import subprocess
import requests
import io
import argparse

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_scores", help="path to input text file with metallo scores")
    parser.add_argument("--input_pdb_rna", help="path to input pdb [RNA] file")    
    parser.add_argument("--input_pdb_mg", help="path to input pdb [Metallo] file")    
    parser.add_argument("--score_cutoff", help="classification score cutoff used when adding MG2+", default = 0.5)    
    parser.add_argument("--output_pdb_prefix", help="path prefix for outputs")
    a = parser.parse_args()
    
    from pymol import cmd
    import pandas as pd
        
    # load train and test data
    scores = pd.read_csv(a.input_scores,delim_whitespace=True, header=None)
    print(scores.head())
    
    # read in pdb
    cmd.load(a.input_pdb_rna, "rna")
    cmd.load(a.input_pdb_mg, "mg")
    
    # loop over scores and set bfactor
    cmd.alter("all", "q = 0.0 ")
    cmd.alter("all", "b = 0.0 ")
    
    for score in scores.itertuples():
        print (score[0], score[1], score[2], score[4])
        cmd.select("pseudo", "mg and index %s and name %s"%(score[1], score[2]))
        cmd.alter("pseudo", "b = '%s'"%(score[4]))
        cmd.alter("pseudo", "q = '%s'"%(score[4]))
        cmd.alter("pseudo", "chain = 'Z'")
        cmd.alter("mg and index %s and name %s"%(score[1], score[2]), "resi = '%s'"%(score[1]))
        
    # remove MG ions that are too close the RNA
    cmd.remove("br. mg within 2 of rna")
    myspace = {'bfact': []}
    cmd.iterate("mg", "bfact.append(b)", space=myspace)    
    bfactors = myspace['bfact']    
    max_score = numpy.max(bfactors)
    
    for score in scores.itertuples():
        print (score[0], score[1], score[2], score[4])
        cmd.select("pseudo", "mg and index %s and name %s"%(score[1], score[2]))
        cmd.alter("pseudo", "b = '%s'"%(score[4]/max_score))
        cmd.alter("pseudo", "q = '%s'"%(score[4]))
        cmd.alter("pseudo", "chain = 'Z'")
        cmd.alter("mg and index %s and name %s"%(score[1], score[2]), "resi = '%s'"%(score[1]))
        
    # save coordinate files    
    cmd.select("best", "br. b>%s and mg"%(a.score_cutoff))    
    cmd.create("complex", "rna mg")
    cmd.save(a.output_pdb_prefix+"_all.pdb", "complex")
    cmd.delete("complex")
    cmd.alter("best", "resn = 'MG'")
    cmd.alter("best", "name = 'MG'")
    cmd.create("complex", "rna best")
    cmd.save(a.output_pdb_prefix+"_best.pdb", "complex")

