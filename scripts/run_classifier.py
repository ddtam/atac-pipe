#!/usr/bin/env python

import sys
import joblib
import pandas as pd

from sklearn.pipeline import Pipeline
from lightgbm.sklearn import LGBMClassifier

pipe_lgbmc_best = joblib.load('src/lgbmc_model.pkl')

data = pd.read_csv(
    sys.argv[1], sep="\t")

celltype_labs = [
    "CLP", "CMP", "GMP", "HSC/MPP", "LMPP", "MEP", "mono", "pDC"
]


def build_proba_matrix(pipe, df):
    output = pd.concat(
        [pd.DataFrame(pipe.predict_proba(df)), pd.DataFrame(pipe.predict(df))],
        axis=1,
        ignore_index=True,
    )

    output.columns = celltype_labs + ["prediction"]

    return output


classification = build_proba_matrix(pipe_lgbmc_best, data)

classification.to_csv(sys.argv[2], sep="\t")
