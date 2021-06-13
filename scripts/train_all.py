import sys
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score
from sklearn.externals import joblib

# Load data
rnas = ["1vt2.txt", "1vvl.txt", "1vvm.txt", "1vvu.txt", "1vvz.txt", "1vw0.txt", "1vx8.txt", "1vx9.txt", "1vy0.txt", "1vy3.txt", "2avy.txt", "2aw4.txt", "2aw7.txt", "2awb.txt", "2i2p.txt", "2i2t.txt", "2i2u.txt", "2i2v.txt", "2j02.txt", "2qal.txt", "2qam.txt", "2qan.txt", "2qao.txt", "2qb9.txt", "2qba.txt", "2qbb.txt", "2qbc.txt", "2qbd.txt", "2qbe.txt", "2qbf.txt", "2qbg.txt", "2qbh.txt", "2qbi.txt", "2qbj.txt", "2qbk.txt", "2qou.txt", "2qov.txt", "2qow.txt", "2qox.txt", "2qoy.txt", "2qoz.txt", "2qp0.txt", "2qp1.txt", "2v46.txt", "2wdg.txt", "2wdh.txt", "2wdk.txt", "2z4k.txt", "2z4l.txt", "2z4m.txt", "2z4n.txt", "3df1.txt", "3df2.txt", "3df3.txt", "3df4.txt", "3i1m.txt", "3i1n.txt", "3i1o.txt", "3i1p.txt", "3i1q.txt", "3i1r.txt", "3i1s.txt", "3i1t.txt", "3i1z.txt", "3i20.txt", "3i21.txt", "3i22.txt", "3i8g.txt", "3i9b.txt", "3i9d.txt", "3kni.txt", "3oaq.txt", "3oar.txt", "3oas.txt", "3oat.txt", "3ofa.txt", "3ofb.txt", "3ofc.txt", "3ofd.txt", "3ofo.txt", "3ofp.txt", "3ofq.txt", "3ofr.txt", "3ofx.txt", "3ofy.txt", "3ofz.txt", "3og0.txt", "3or9.txt", "3ora.txt", "3orb.txt", "3r8s.txt", "3r8t.txt", "3v22.txt", "3v24.txt", "3v25.txt", "3v26.txt", "3v27.txt", "3v28.txt", "3v29.txt", "3v2c.txt", "3v2d.txt", "3v2e.txt", "3v2f.txt", "3v6u.txt", "4a17.txt", "4a1a.txt", "4a1c.txt", "4a1e.txt", "4bpe.txt", "4bpn.txt", "4bpo.txt", "4bpp.txt", "4bye.txt", "4cr1.txt", "4dh9.txt", "4dha.txt", "4dhb.txt", "4dhc.txt", "4g5m.txt", "4g5t.txt", "4g5v.txt", "4gaq.txt", "4gar.txt", "4gas.txt", "4gau.txt", "4gd1.txt", "4gd2.txt", "4hub.txt", "4kix.txt", "4kiy.txt", "4kiz.txt", "4kj0.txt", "4kj1.txt", "4kj2.txt", "4kj3.txt", "4kj4.txt", "4kj5.txt", "4kj6.txt", "4kj7.txt", "4kj8.txt", "4kj9.txt", "4kja.txt", "4kjb.txt", "4kjc.txt", "4kx2.txt", "4nvu.txt", "4nvv.txt", "4nvw.txt", "4nvx.txt", "4nvy.txt", "4nvz.txt", "4nw0.txt", "4nw1.txt", "4pe9.txt", "4pea.txt", "4peb.txt", "4pec.txt", "4qcm.txt", "4qcn.txt", "4qco.txt", "4qcp.txt", "4qcq.txt", "4qcr.txt", "4qcs.txt", "4qct.txt", "4qcu.txt", "4qcv.txt", "4qcw.txt", "4qcx.txt", "4qcy.txt", "4qcz.txt", "4qd0.txt", "4qd1.txt", "4tol.txt", "4tom.txt", "4ton.txt", "4too.txt", "4tou.txt", "4tov.txt", "4tow.txt", "4tox.txt", "4tp0.txt", "4tp1.txt", "4tp2.txt", "4tp3.txt", "4tp4.txt", "4tp5.txt", "4tp6.txt", "4tp7.txt", "4tp8.txt", "4tp9.txt", "4tpa.txt", "4tpb.txt", "4tpc.txt", "4tpd.txt", "4tpe.txt", "4tpf.txt", "3i56.txt", "3g71.txt", "1hq1.txt", "3ud3.txt", "3ccl.txt", "2b8s.txt", "1yij.txt", "3mxh.txt", "3ccm.txt", "1k8a.txt", "4ji8.txt", "3oxj.txt", "1dfu.txt", "3gx5.txt", "4p5j.txt", "3po3.txt", "1zbl.txt", "3mum.txt", "2nug.txt", "2nuf.txt", "1nji.txt", "3po2.txt", "3g6e.txt", "4enb.txt", "1zbi.txt", "4m4o.txt", "3oxm.txt", "3ccj.txt", "3d2x.txt", "4a93.txt", "1yjw.txt", "1w2b.txt", "2uuc.txt", "4b3r.txt", "3ud4.txt", "3u2e.txt", "4en5.txt", "4lfa.txt", "4enc.txt", "4lf7.txt", "4ena.txt", "4lfc.txt", "3dil.txt", "2uua.txt", "3aoh.txt", "4tzx.txt", "3cxc.txt", "2r1s.txt", "3oww.txt", "4jf2.txt", "2oiy.txt", "1m90.txt", "2qbz.txt", "4lfb.txt", "2hw8.txt", "1jzv.txt", "4lx5.txt", "1xjr.txt", "1vq6.txt", "3ucu.txt", "1kd1.txt", "2g91.txt", "1nuj.txt", "3hvr.txt", "3td0.txt", "4dv7.txt", "4duz.txt", "3ho1.txt", "4fb0.txt", "4a3c.txt", "3hm9.txt", "4a3b.txt", "1y99.txt", "1yls.txt", "4ohz.txt", "4dv6.txt", "1ytu.txt", "2yie.txt", "2qex.txt", "3rer.txt", "1vq7.txt", "1m1k.txt", "354d.txt", "2c4r.txt", "1q86.txt", "4bxx.txt", "2h0s.txt", "1vq5.txt", "2a43.txt", "2ez6.txt", "3cgs.txt", "4dv4.txt", "4tra.txt", "3hax.txt", "4duy.txt", "1pjo.txt", "1n8r.txt", "1kqs.txt", "4dv5.txt", "3htx.txt", "2gcv.txt", "2yif.txt", "2g8v.txt", "437d.txt", "1vq4.txt", "1q7y.txt", "1ehz.txt", "4dv1.txt", "1n78.txt", "4a3e.txt", "4a3d.txt", "4dv0.txt", "2gcs.txt", "4jnx.txt", "4m2z.txt", "3r1l.txt", "2h0w.txt", "1vqp.txt", "1q82.txt", "4ifd.txt", "1fuf.txt", "1qvg.txt", "2ao5.txt", "2g8f.txt", "1lng.txt", "4dv2.txt", "4a3f.txt", "4dv3.txt", "1l8v.txt", "4ari.txt", "2q1r.txt", "2o5i.txt", "1qvf.txt", "1dno.txt", "1q81.txt", "2z75.txt", "3iin.txt", "2otj.txt", "4un3.txt", "3bo3.txt", "3mei.txt", "1ik5.txt", "4a3k.txt", "2r8s.txt", "3zgz.txt", "4m30.txt", "1dnt.txt", "2h0x.txt", "3bo2.txt", "2fqn.txt", "1qu2.txt", "1d4r.txt", "1s72.txt", "1tez.txt", "3jxq.txt", "462d.txt", "4faw.txt", "4hor.txt", "3cd6.txt", "1nuv.txt", "2g8h.txt", "1vqk.txt", "3la5.txt", "4pqv.txt", "1tra.txt", "2otl.txt", "4by1.txt", "1vqo.txt", "4nxn.txt", "4un5.txt", "1vq9.txt", "3b4c.txt", "3ucz.txt", "2g5k.txt", "1egk.txt", "4a3l.txt", "4far.txt", "4a3m.txt", "4qln.txt", "3cpw.txt", "3b4b.txt", "1vq8.txt", "4un4.txt", "1vqn.txt", "4bwm.txt", "1vql.txt", "4nxm.txt", "3egz.txt", "4j7l.txt", "1ntb.txt", "4oji.txt", "4qlm.txt", "4q9q.txt", "3b4a.txt", "2zxu.txt", "1vqm.txt", "1o3z.txt", "1dul.txt", "1x8w.txt", "3p49.txt", "4lf9.txt", "1mms.txt", "3oxb.txt", "3cce.txt", "1k9m.txt", "3ccr.txt", "4ji0.txt", "3owz.txt", "4dr5.txt", "2nz4.txt", "1yj9.txt", "3v7e.txt", "1yhq.txt", "3hk2.txt", "1yjn.txt", "4dr4.txt", "4ji1.txt", "3ccs.txt", "3cc2.txt", "3d2v.txt", "3zc0.txt", "4lf8.txt", "3nkb.txt", "4kzd.txt", "3wfr.txt", "1ffk.txt", "2oe5.txt", "3l0u.txt", "1l2x.txt", "1k73.txt", "3s14.txt", "2bx2.txt", "3dd2.txt", "4ji3.txt", "3ccq.txt", "4dr6.txt", "1y26.txt", "2r20.txt", "3mur.txt", "4dr7.txt", "4ji2.txt", "3s15.txt", "1jj2.txt", "6tna.txt", "1feu.txt", "1f27.txt", "1ffy.txt", "1hr2.txt", "3oxd.txt", "4ji6.txt", "1evv.txt", "4dr3.txt", "3q3z.txt", "3q1r.txt", "3muv.txt", "1yi2.txt", "3hhn.txt", "3hjf.txt", "4dr2.txt", "2ppb.txt", "3ccu.txt", "3cc4.txt", "4ji7.txt", "3d2g.txt", "1kc8.txt", "3oxe.txt", "4e8m.txt", "2ho6.txt", "4ji5.txt", "3twh.txt", "3ssf.txt", "3q1q.txt", "4oog.txt", "3mut.txt", "4dr1.txt", "3owi.txt", "2gdi.txt", "4ji4.txt", "3ccv.txt", "3cc7.txt", "3ox0.txt", "3ivk.txt", "4e8n.txt", "1j1u.txt", "2c0b.txt", "4p9r.txt","4g6r.txt", "3v23.txt", "1vs5.txt", "1vs6.txt", "1vs7.txt", "1vs8.txt"]

# For each RNA, put it in the test set, along with the RNAs that have >80% sequence similarity
for rna in rnas:
	training_set_files = []
	test_set_files = []
	pdbID = rna[0:4]
	print(pdbID)
	#os.chdir("/home/zhang/Desktop/Mg_project/data/rna_sequences/" + pdbID)
	test_set_files.append(rna)
	f = open("../data/rna_sequences/" + pdbID + "/" + pdbID + "_list.txt")
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
		    rna_file = pd.read_csv("../results_s_2.8/"+training_set_rnas, delim_whitespace=True, engine='python')
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
		    rna_file = pd.read_csv("../results_s_2.8/"+test_set_rnas, delim_whitespace=True, engine='python')
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
	if os.path.exists("../data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt"):
		os.remove("../data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt")
	f = open("../data/rna_sequences/" + pdbID + "/" + pdbID + "_report.txt", "a")
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
	path = "../data/rna_sequences/" + pdbID + "/" + pdbID + "_predictor.pkl"
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
