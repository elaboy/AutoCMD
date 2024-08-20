from typing import List
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from Objects import engine, session, PSM
from pipe import *
from utils import *

def charge_state_distribution(charge_states:List[int]) -> None:
    fig, ax = plt.subplots()

    ax.hist(charge_states)
    # ax.title = "Precursor Charge Distribution"

def pairwise_charge_state_distribution(charge_states_one:List[int],
                                       charge_states_two:List[int]) -> None:
    fig, ax = plt.subplots()

    ax.hist(charge_states_one)
    ax.hist(charge_states_two)

def HI_comparer(dataframes: List[pd.DataFrame], labels:List[str]) -> None:
    fig, ax = plt.subplots(len(dataframes))

    for idx, dataframe in enumerate(dataframes):
        dataframe = dataframe.sort_values(by=" Scan Retention Time ", ascending=True)
        
        dataframe[" Scan Retention Time "] = (dataframe[" Scan Retention Time "]-
                                              dataframe[" Scan Retention Time "].mean()) / \
                                                dataframe[" Scan Retention Time "].std()
        dataframe[" ChronologerHI"] = (dataframe[" ChronologerHI"]-
                                        dataframe[" ChronologerHI"].mean()) / \
                                        dataframe[" ChronologerHI"].std()
        
        ax[idx].scatter(dataframe[" Scan Retention Time "],
                        dataframe[" ChronologerHI"],
                        label = labels[idx],
                        s = 0.3)
        
    plt.scatter(len(dataframes[0]) - 1, dataframes[0][" Scan Retention Time "].values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize your results")

    parser.add_argument("--data",
                        metavar='data',
                        nargs="*",
                        type=str,
                        help="list paths to the data")
    
    parser.add_argument("--labels",
                        metavar="labels",
                        nargs="*",
                        type=str,
                        help="list the labels for each dataset in the same order")
    
    # parse the arguments
    args = parser.parse_args()
    
    # save as variables
    directories = args.data
    labels = args.labels

    #inputs for analysis
    analysis_to_do = input("""
    Type the analysis to be made (the integer) and press the Enter key: \n
        0 -> all analysis
        1 -> compare Chronologer HI and Scan Reported Retention Time
    """)

    # TODO: make a method for all analysis to run them in a Pool

    match analysis_to_do:
        case "1":
            print("crunching numbers for you :) ...")
            
            dataframes = []

            for path in directories:
                dataframes.append(pd.read_csv(path, sep="\t", low_memory=False))

            HI_comparer(dataframes=dataframes, labels=labels)

    plt.show()