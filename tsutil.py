import numpy as np
import pandas as pd


enum_times = [
 '3wk',
 '4wkTC',
 '4wkAB',
 '6wkTC',
 '6wkAB',
 '8wkTC',
 '8wkAB',
 '10wkTC',
 '10wkAB']


enum_times_TC = [
 '3wk',
 '4wkTC',
 '6wkTC',
 '8wkTC',
 '10wkTC']


enum_times_AB = [
 '3wk',
 '4wkAB',
 '6wkAB',
 '8wkAB',
 '10wkAB']


def read_ts_data(CV_limit: float):
    mb = pd.read_hdf("data/generated/biom_ts.hdf5", "CV" + str(CV_limit).replace('.', ''))
    label = pd.read_hdf("data/generated/label_ts.hd5", "default")
    print("Using {} metabolome and {} bacteria data".format(len(mb), len(label)))

    mb_TC = mb.filter(enum_times_TC)
    mb_AB = mb.filter(enum_times_AB)
    label_TC = label.filter(enum_times_TC)
    label_AB = label.filter(enum_times_AB)

    mb_TC_np = mb_TC.to_numpy().astype(float)
    mb_AB_np = mb_AB.to_numpy().astype(float)
    label_TC_np = label_TC.to_numpy().astype(float)
    label_AB_np = label_AB.to_numpy().astype(float)

    mb_np = mb.to_numpy()[:, 1:].astype(float)
    label_np = label.to_numpy()[:, 1:].astype(float)

    return mb_np, label_np, mb_TC_np, label_TC_np, mb_AB_np, label_AB_np


def save_score(score: str, type: str, reg: str):
    with open(f"scores_{reg}_{type}", 'a') as f:
        f.write(str(score))






    
    