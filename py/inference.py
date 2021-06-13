import joblib
import sys
import os
from train_utils import inference

def main():
    path = sys.argv[1]
    home = sys.argv[2]
    feathome, rna = os.path.dirname(path), os.path.basename(path).split(".")[0]
    # set path
    modelhome = "%s/model" %home
    outhome = "predict_output"

    model = "%s/s28_model_no_cut.pkl" %modelhome
    clf = joblib.load(model)

    inference(clf, [rna], feathome, outhome, drop_M=False)

    return 0

if __name__ == '__main__':
    main()
