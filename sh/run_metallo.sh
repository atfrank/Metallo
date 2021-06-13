#!/bin/bash

if [[ $# -ne 2 ]]
then
	echo "usage: $0 <path-to-pdb-file> <path-to-metallo-repo>"
else

	# initialize variables
	pdb_file=$1
	metallo=$2
	rna="user"

	# clean up old results
	rm -rfv output; mkdir output
	cd output
	mkdir results addatom predict_output

	# copy need files
	ln -s ${metallo}/py/train_utils.py .

	# make a fix names and make local copy
	sed 's/GUA/  G/g' ../${pdb_file} | \
	sed 's/ADE/  A/g' | \
	sed 's/CYT/  C/g' | \
	sed 's/URA/  U/g' | \
	sed 's/RG3/  G/g' | \
	sed 's/RA3/  A/g' | \
	sed 's/RC3/  C/g' | \
	sed 's/RU3/  U/g' | \
	sed 's/RG5/  G/g' | \
	sed 's/RA5/  A/g' | \
	sed 's/RC5/  C/g' | \
	sed 's/RU5/  U/g' | \
	sed 's/RG/ G/g' | \
	sed 's/RA/ A/g' | \
	sed 's/RC/ C/g' | \
	sed 's/RU/ U/g' | \
	sed 's/O1P/OP1/g' | \
	sed 's/O2P/OP2/g' | tee ${rna}.pdb

	# featurize
	echo "featurizing..."	
	s=2.8 # spacing
	sel="O -sel OP1 -sel OP2 -sel O2' -sel O3' -sel O4' -sel O5' -sel N1 -sel N2 -sel N3 -sel N4 -sel N6 -sel N7 -sel N9 -sel O6 -sel O2 -sel O4"
	${metallo}/bin/metallo -dmyf -s $s -sel $sel  -row J -row L -row M ${rna}.pdb
	
	# infer
	echo "infering..."
	python ${metallo}/py/inference.py results/${rna}.txt ${metallo}

	# get list of features
	paste <(awk '{if(NR!=1) print $1, $2, $3}' results/${rna}.txt) <(cat predict_output/${rna}.txt) > predict_output/pred.txt

	# create pdb files
	python ${metallo}/py/addtoPDB.py --input_scores predict_output/pred.txt --input_pdb_rna ${rna}.pdb --input_pdb_mg addatom/addatom_${rna}.pdb --output_pdb addatom/scores_addatom_${rna} --score_cutoff 0.95	

	# add citation info
	touch CITATION

	# remove file
	rm -rfv train_utils.py predict_output results
	mv addatom metallo-pdbs

	echo "Metallo: A Machine Learning Tool for Identifying Magnesium Binding Sites in RNA" >> CITATION
	echo "Jingru Xie, Lichirui Zhang, and Aaron T. Frank, 2020 (in preparation)" >> CITATION
	echo "Source code: https://github.com/atfrank/Metallo" >> CITATION

fi