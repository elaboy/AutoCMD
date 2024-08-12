from typing import List
from Objects import PSM
import matplotlib.pyplot as plt

def charge_state_distributions(data:List[List[PSM]], labels:List[str]) -> None:
    fig, ax = plt.subplots()

    for index in len(data):
        charge_states = []
        for i in data[index]:
            charge_states.append(i.precursor_charge)
        ax.hist(charge_states, label=labels[index])
    
    plt.legend()
    return None