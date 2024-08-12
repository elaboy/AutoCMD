from typing import List
import matplotlib.container
import matplotlib.figure
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from enum import Enum
import os
import sqlalchemy as sa
from sqlalchemy.orm import Mapped, mapped_column, sessionmaker, declarative_base
import argparse
from Objects import database, session, PSM


def charge_state_distribution(charge_states:List[int]) -> None:
    fig, ax = plt.subplots()

    ax.hist(charge_states)
    # ax.title = "Precursor Charge Distribution"

def pairwise_charge_state_distribution(charge_states_one:List[int], charge_states_two:List[int]) -> None:
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
    
    parser.add_argument("--in_memory",
                        metavar="in_memory",
                        type=str,
                        help="Indicateds if the database will be in memory or it should be saved. Y for in memory, N for saving")
    
    # parse the arguments
    args = parser.parse_args()
    
    # save as variables
    directories = parser.search_task_directories
    labels = parser.labels

    df = pd.read_csv(os.path.join(directory, "AllPeptides.psmtsv"), sep="\t")
    df = df.rename_axis("id").reset_index()
    df.to_sql(con=database, name=PSM.__tablename__, if_exists="append", index=True)
    
    results = session.query(PSM).all()
    charge_state_distribution([i.precursor_charge for i in results])
    charge_state_distribution([i.precursor_charge for i in results])
    
    plt.show()