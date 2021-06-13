# Metallo
Metallo: A Machine Learning Tool for Classifying Magnesium Binding Sites in RNA


## Quick Start
```
git clone https://github.com/atfrank/Metallo.git
cd Metallo/
make clean
make
```
### Install Dependencies

#### Python Modules
```
conda create -n metallo
conda activate metallo

conda install -c schrodinger pymol pandas -y 
pip install scikit-learn==0.22.2.post1 joblib matplotlib
```

### Using Metallo
* Main script is `./sh/run_metallo.sh`. 
```
usage: ./sh/run_metallo.sh <path-to-pdb-file> <path-to-metallo-repo>
```

* Example:
```
conda activate metallo
./sh/run_metallo.sh test/SAM.pdb  ~/Documents/GitHub/Metallo
```
This generates a coordinate (PDB formatted) files stored in `output/metallo-pdbs/`.


* `output/metallo-pdbs/scores_addatom_user_best.pdb` : original coordinate file
* `output/metallo-pdbs/addatom_user.pdb`: original + grid points coordinate file 
* `output/metallo-pdbs/scores_addatom_user_all.pdb`: original + highest scoring grid points coordinate file

## Metallo via SMALTR
* Metallo can also be accessed via: [SMALTR](http://smaltr.org/).


## Publications
* Jingru Xie, Lichirui Zhang, and Aaron T. Frank. "Metallo: A Machine Learning Tool for Identifying Magnesium Binding Sites in RNA" (in preparation)


## COMMERCIAL USE LICENSE: 

If you are interested in commercial licensing of these applications (clinical, operational, etc.) please contact the University of Michigan Office of Technology Transfer for a quote and licensing options.

Drew Bennett - https://techtransfer.umich.edu/team/drew-bennett/

or

techtransfer@umich.edu





