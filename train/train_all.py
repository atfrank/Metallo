import sys
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score
from sklearn.externals import joblib

# Load data
rnas = ["1ehz.txt", "1evv.txt", "1feu.txt", "1ffy.txt", "1fuf.txt", "1hr2.txt", "1jj2.txt", "1k73.txt", "1k9m.txt", "1kd1.txt", "1kqs.txt", "1l2x.txt", "1lng.txt", "1m1k.txt", "1m90.txt", "1mms.txt", "1n78.txt", "1nji.txt", "1ntb.txt", "1nuj.txt", "1nuv.txt", "1pjo.txt", "1q7y.txt", "1q81.txt", "1q82.txt", "1qu2.txt", "1qvf.txt", "1qvg.txt", "1s72.txt", "1tra.txt", "1vq4.txt", "1vq5.txt", "1vq6.txt", "1vq7.txt", "1vq8.txt", "1vq9.txt", "1vqk.txt", "1vql.txt", "1vqo.txt", "1vqp.txt", "1x8w.txt", "1xjr.txt", "1yhq.txt", "1yi2.txt", "1yij.txt", "1yj9.txt", "1yjn.txt", "1yls.txt", "2a43.txt", "2gcs.txt", "2gcv.txt", "2gdi.txt", "2h0s.txt", "2h0w.txt", "2h0x.txt", "2ho6.txt", "2nz4.txt", "2otj.txt", "2otl.txt", "2qbz.txt", "2qex.txt", "2r8s.txt", "2yie.txt", "2yif.txt", "2z75.txt", "2zxu.txt", "354d.txt", "3b4a.txt", "3b4b.txt", "3b4c.txt", "3bo3.txt", "3cc2.txt", "3cc7.txt", "3ccj.txt", "3ccl.txt", "3ccm.txt", "3ccr.txt", "3ccu.txt", "3ccv.txt", "3cxc.txt", "3d2g.txt", "3d2v.txt", "3d2x.txt", "3dd2.txt", "3dil.txt", "3egz.txt", "3g71.txt", "3gx5.txt", "3hhn.txt", "3ivk.txt", "3l0u.txt", "3la5.txt", "3mum.txt", "3mut.txt", "3mxh.txt", "3nkb.txt", "3owi.txt", "3oww.txt", "3owz.txt", "3ox0.txt", "3oxb.txt", "3oxd.txt", "3oxm.txt", "3p49.txt", "3rer.txt", "3u2e.txt", "3ucz.txt", "3v7e.txt", "3zgz.txt", "462d.txt", "4dr2.txt", "4dr3.txt", "4dr4.txt", "4dr5.txt", "4dr6.txt", "4dr7.txt", "4duy.txt", "4duz.txt", "4dv0.txt", "4dv1.txt", "4dv2.txt", "4dv4.txt", "4dv7.txt", "4e8m.txt", "4e8n.txt", "4ena.txt", "4enb.txt", "4enc.txt", "4faw.txt", "4ifd.txt", "4jf2.txt", "4ji0.txt", "4ji1.txt", "4ji2.txt", "4ji3.txt", "4ji4.txt", "4ji5.txt", "4ji6.txt", "4ji7.txt", "4ji8.txt", "4kzd.txt", "4lf7.txt", "4lf8.txt", "4lf9.txt", "4lfb.txt", "4m30.txt", "4m4o.txt", "4nxm.txt", "4nxn.txt", "4oji.txt", "4pqv.txt", "4q9q.txt", "4qlm.txt", "4tra.txt", "4tzx.txt", "6tna.txt"]

# For each RNA, put it in the test set, along with the RNAs that have >80% sequence similarity
for rna in rnas:
	training_set_files = []
	test_set_files = []
	pdbID = rna[0:4]
	print(pdbID)
	test_set_files.append(rna)
	f = open("data/rna_sequences/" + pdbID + "/" + pdbID + "_list.txt")
	for x in f:
		x = x[0:4] + ".txt"
		test_set_files.append(x.lower())
	f.close()
	
	for x in rnas:
		if x not in (name.lower() for name in test_set_files):
			training_set_files.append(x.lower())
	X = []
	Y = []
	Y_test = []
	X_test = []

	for training_set_rnas in training_set_files:
		try:
		    rna_file = pd.read_csv("results_s_2.8/"+training_set_rnas, delim_whitespace=True, engine='python')
		except IOError:
			print(training_set_rnas)
			continue
		rna_rawdata = pd.DataFrame(rna_file)
		row = rna_rawdata.shape[0]
		col = rna_rawdata.shape[1]
		for i in range(0,row):
		    if rna_rawdata.iloc[i,1] == "J":
		        Y.append(0)
		        X.append(np.asarray(rna_rawdata.iloc[i,3:col-1]))
		    if rna_rawdata.iloc[i,1] == "L":
		        Y.append(1)
		        X.append(np.asarray(rna_rawdata.iloc[i,3:col-1]))

	for test_set_rnas in test_set_files:
		try:
		    rna_file = pd.read_csv("results_s_2.8/"+test_set_rnas, delim_whitespace=True, engine='python')
		except IOError:
			print(test_set_rnas)
			continue
		rna_rawdata = pd.DataFrame(rna_file)
		row = rna_rawdata.shape[0]
		for i in range(0,row):
		    if rna_rawdata.iloc[i,1] == "J":
		        Y_test.append(0)
		        X_test.append(np.asarray(rna_rawdata.iloc[i,3:col-1]))
		    if rna_rawdata.iloc[i,1] == "L":
		        Y_test.append(1)
		        X_test.append(np.asarray(rna_rawdata.iloc[i,3:col-1]))


	X = np.asarray(X)
	Y = np.asarray(Y)
	X_test = np.asarray(X_test)
	Y_test = np.asarray(Y_test)
	print("prepare data done!")

	#build a RF classifier
	clf = RandomForestClassifier(n_estimators=100)
	clf.fit(X, Y)

	predict = clf.predict(X_test)
	print(predict)
	

	#report performance
	if os.path.exists("data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt"):
		os.remove("data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt")
	f = open("data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt", "a")
	f.write("f1_score precision recall test_acc train_acc\n")
	f.write("{:.8f} ".format(f1_score(Y_test, predict)))
	f.write("{:.8f} ".format(precision_score(Y_test, predict)))
	f.write("{:.8f} ".format(recall_score(Y_test, predict)))
	f.write("{:.8f} ".format(clf.score(X_test,Y_test)))
	f.write("{:.8f}\n\n".format(clf.score(X,Y)))
	f.write(np.array2string(confusion_matrix(Y_test, predict), separator=', '))
	f.write("\n")
	for item in test_set_files:
		f.write("%s\n" % item)
	f.close()
	path = "data/rna_sequences/" + pdbID + "/" + pdbID + "_predictor.pkl"
	joblib.dump(clf, path, compress=3)
	print("f1_score {:.8f}".format(f1_score(Y_test, predict)))
	print("precision {:.8f}".format(precision_score(Y_test, predict)))
	print("recall {:.8f}".format(recall_score(Y_test, predict)))
	print("test_acc {:.8f}".format(clf.score(X_test,Y_test)))
	print("train_acc {:.8f}".format(clf.score(X,Y)))
	print("\n")
	print(confusion_matrix(Y_test, predict))
	print("\n")
	
	del predict, X, Y, X_test, Y_test, clf, training_set_files, test_set_files
