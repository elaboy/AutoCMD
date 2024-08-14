from typing import List
from Objects import PSM
import matplotlib.pyplot as plt
import pandas as pd
import chronologer.Predict_RT
from Objects import session

def charge_state_distributions(data, labels:List[str]) -> None:
    fig, ax = plt.subplots()

    for i in data:
        print(i)
        break

    for index in len(data):
        charge_states = []
        for i in data[index]:
            charge_states.append(i.precursor_charge)
        ax.hist(charge_states, label=labels[index])
    
    plt.legend()
    return None

# def get_chronologer_HI(full_sequences:List[str]) -> pd.DataFrame:
#     df = chronologer.Predict_RT.predict_HI(full_sequences)
#     return df 