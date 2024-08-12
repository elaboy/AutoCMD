from typing import List
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
from sqlalchemy.orm import session, Session, sessionmaker
from Objects import *
from enum import Enum

class ResultType(Enum):
    AllPeptides = 1
    AllQuantifiedPeaks = 2
    AllPSMs = 3
    IndividualPeptides = 4
    IndividualQuantifiedPeaks = 5
    IndividualPSMs = 6

def reader(path:str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")

def get_dataframes(paths:List[str]) -> List[pd.DataFrame]:

    if len(paths) > mp.cpu_count() -1:
        pool = Pool(mp.cpu_count()-1)
    else:
        pool = Pool(len(paths))

    df_list = pool.map(reader, paths)
    
    for df in df_list:
        df = df.rename_axis("id").reset_index()

    return df_list

def setup_database(dataframe:pd.DataFrame, result_type:ResultType, if_exists="append", index=True) -> None:
    #TODO: finish method
    
    return None

def setup_databases(dataframes:List[pd.DataFrame], results_type:List[ResultType], if_exists:List[str], index:List[Boolean]) -> None:
    #TODO: finish method
    
    return None

def get_filtered_data(session:Session, table_name:str, ambiguity_level:str, q_value_threshold:float,
              pep_threshold:float, pep_q_value_threshold:float) -> List[PSM]:
    session.query()


