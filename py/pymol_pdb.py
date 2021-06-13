import os
import argparse

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_pdb", help="input pdb path")    
    parser.add_argument("--output_pdb", help="output pdb path")    
    a = parser.parse_args()
    
    from pymol import cmd    
    cmd.load(a.input_pdb, "rna")

    # clean up
    cmd.remove("rna and resn RG3 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RA3 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RC3 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RU3 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RG5 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RA5 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RC5 and name P+OP1+OP2+OP3")
    cmd.remove("rna and resn RU5 and name P+OP1+OP2+OP3")
    cmd.remove("hydrogen")
    
    cmd.alter("resn rC+C+RC+RC3+RC5", "resn = 'CYT'")
    cmd.alter("resn rA+A+RA+RA3+RA5", "resn = 'ADE'")
    cmd.alter("resn rU+U+RU+RU3+RU5", "resn = 'URA'")
    cmd.alter("resn rG+G+RG+RG3+RG5", "resn = 'GUA'")    
    cmd.save(a.output_pdb, "rna")

