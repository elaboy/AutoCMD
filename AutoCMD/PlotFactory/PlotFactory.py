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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize your results")

    parser.add_argument("--search_task_directories",
                        metavar='search_task_directories',
                        nargs="*",
                        type=str,
                        help="list search task directories to the data")
    
    parser.add_argument("--labels",
                        metavar="labels",
                        nargs="*",
                        type=str,
                        help="list the labels for each dataset in the same order")
    
    # parse the arguments
    args = parser.parse_args()
    
    # save as variables
    directories = args.search_task_directories
    labels = args.labels

    #inputs for analysis
    analysis_to_do = input("""
    Type the analysis to be made (the integer) and press the Enter key: \n
        0 -> all analysis
        1 -> charge state distributions
    """)

    match analysis_to_do:
        case "1":
            print("crunching numbers for you :) ...")
            
            dataframes = get_dataframes(directories)
            
            setup_database(dataframe=dataframes[0], 
                            result_type=ResultType.AllPeptides, 
                            if_exists="append",
                            index=True)
            
            print(len(session.query(PSM).all()))
            # results = get_PSM_filtered_data(ambiguity_level="1", 
            #                                 q_value_threshold=0.01, 
            #                                 pep_threshold=0.5, 
            #                                 pep_q_value_threshold=0.01)
            
            # charge_state_distributions(data = results, labels = labels)
    plt.show()