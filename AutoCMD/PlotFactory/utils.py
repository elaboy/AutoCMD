from typing import List
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
from Objects import *
from enum import Enum
from sqlalchemy import select

class ResultType(Enum):
    AllPeptides = 1
    AllQuantifiedPeaks = 2
    AllPSMs = 3
    IndividualPeptides = 4
    IndividualQuantifiedPeaks = 5
    IndividualPSMs = 6

def get_result_type(result_type: ResultType) -> str:
    match result_type:
        
        case ResultType.AllPeptides | ResultType.AllPSMs:
            return PSM.__tablename__
        
        case ResultType.AllQuantifiedPeaks:
            return QuantifiedPeak.__tablename__
        
    return None

def reader(path:str, label:str) -> pd.DataFrame:
    df =  pd.read_csv(path, sep="\t", low_memory=False)
    
    if "Precursor Intensity" in df:
        print("Dataframe check: Pass")
    else:
        print("Dataframe check: Failed, looks like an old MM version. Fixing...")
        df["Precursor Intensity"] = 0
    
    df["Label"] = label
    
    return df

#TODO: address mismatch in column size
def get_dataframes(paths:List[str], labels:List[str]) -> List[pd.DataFrame]:

    if len(paths) > mp.cpu_count() -1:
        pool = Pool(mp.cpu_count()-1)
    else:
        pool = Pool(len(paths))

    df_list = pool.starmap(reader, zip(paths, labels))
    
    for df in df_list:
        df = df.rename_axis("id")

    return df_list

def setup_database(dataframe:pd.DataFrame, result_type:ResultType, if_exists="append", index=True) -> None:
    print(dataframe.head())
    dataframe.to_sql(name="psms",con=engine, if_exists=if_exists, index=index, index_label="id")
    return None

def setup_databases(dataframes:List[pd.DataFrame], results_type:List[ResultType],
                     if_exists:List[str], index:List[Boolean]) -> None:
    for idx, df in enumerate(dataframes):
        setup_database(df, results_type[idx], if_exists[idx], index[idx])

    return None

def get_PSM_filtered_data(ambiguity_level:str, q_value_threshold:float,
              pep_threshold:float, pep_q_value_threshold:float) -> List[PSM]:
    
    results = session.query(PSM).filter(PSM.ambiguity_level == ambiguity_level, 
                                         PSM.q_value <= q_value_threshold,
                                         PSM.posterior_error_probability <= pep_threshold,
                                         PSM.posterior_error_probability_q_value <= pep_q_value_threshold,
                                         PSM.decoy == "N",
                                         PSM.decoy_contaminant_target == "T").all()
    
    return results

#TODO: method to predict HI from full sequences using chronologer

