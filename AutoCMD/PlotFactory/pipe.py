from typing import List
from Objects import PSM
import matplotlib.pyplot as plt
import pandas as pd
import chronologer.Predict_RT
from Objects import session
import itertools

def charge_state_distributions(results:List[PSM], labels:List[str]) -> None:
    fig, ax = plt.subplots()

    data = []
    for label in labels:
        data.append([int(psm.precursor_charge) for psm in results if psm.label == label])

    for idx, label in enumerate(labels):
        ax.hist(data[idx], label = label)

    plt.legend()
    return None