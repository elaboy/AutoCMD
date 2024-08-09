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
from Objects import PSM, QuantifiedPeak

database = sa.create_engine('sqlite:///:memory:')
Session = sessionmaker(bind=database)
session = Session()
Base = declarative_base()


class PolymerType(Enum):
    Peptide = 1
    DNA = 2
    RNA = 3

class Polymer:
    def __init__(self, polymer_name:str, polymer_type:PolymerType, polymer_sequence:str) -> None:
        self.polymer_name = polymer_name
        self.polymer_type = polymer_type
        self.polymer_sequence = polymer_sequence

    def __str__(self) -> str:
        return f"Polymer Name: {self.polymer_name}, Polymer Type: {self.polymer_type}, Polymer Sequence: {self.polymer_sequence}"

class QuantifiedPeptide:
    def __init__(self, row:pd.Series):
        self.sequence = row.get("Sequence")
        self.base_sequence = row.get("Base Sequence")
        self.protein_groups = row.get("Protein Groups")
        self.gene_names = row.get("Gene Names")
        self.organism = row.get("Organism")
        self.file_intensity = row[5]
        self.detection_type = row[6]

class QuantifiedProteinGroup:
    def __init__(self, row:pd.Series):
        self.protein_accession = row.get("Protein Accession")
        self.gene = row.get("Gene")
        self.organism = row.get("Organism")
        self.protein_full_name = row.get("Protein Full Name")
        self.protein_unmodified_mass = row.get("Protein Unmodified Mass")
        self.number_of_peptides_in_group = row.get("Number of Peptides in Group")
        self.unique_peptides = row.get("Unique Peptides")
        self.shared_peptides = row.get("Shared Peptides")
        self.number_of_unique_peptides = row.get("Number of Unique Peptides")
        self.sequence_coverage_fraction = row.get("Sequence Coverage Fraction")
        self.sequence_coverage = row.get("Sequence Coverage")
        self.sequence_coverage_with_mods = row.get("Sequence Coverage with Mods")
        self.fragment_sequence_coverage = row.get("Fragment Sequence Coverage")
        self.modification_info_list = row.get("Modification Info List")
        self.file_intensity = row[14]
        self.number_of_psms = row.get("Number of PSMs")
        self.protein_decoy_contaminant_target = row.get("Protein Decoy/Contaminant/Target")
        self.protein_cumulative_target = row.get("Protein Cumulative Target")
        self.protein_cumulative_decoy = row.get("Protein Cumulative Decoy")
        self.protein_q_value = row.get("Protein QValue")
        self.best_peptide_score = row.get("Best Peptide Score")
        self.best_peptide_notch_q_value = row.get("Best Peptide Notch QValue")

class Results:
    def __init__(self, search_task_path:str) -> None:
        self.search_task_path = search_task_path
        self.all_psms_dataframe = self.__read_all_psms_file()
        self.all_peptides_dataframe  = self.__read_all_peptides_file()
        self.all_quantified_peaks_dataframe = self.__read_all_quantified_peaks_file()
        self.results_dataframe = self.__read_results_file()
        self.experimental_design_dataframe = self.__read_experimental_design_file()

    def __read_all_psms_file(self):
        try:
            all_psms_path = self.search_task_path + 'AllPSMs.psmtsv'
            all_psms_dataframe = pd.read_csv(all_psms_path, sep='\t')
            return all_psms_dataframe
        except Exception as e:
            print("Error reading AllPSMs.psmtsv file: ", e)
        
        return None
    
    def __read_all_peptides_file(self):
        try:
            all_peptides_path = self.search_task_path + 'AllPeptides.psmtsv'
            all_peptides_dataframe = pd.read_csv(all_peptides_path, sep='\t')
            return all_peptides_dataframe
        except Exception as e:
            print("Error reading AllPeptides.psmtsv file: ", e)

        return None
    
    def __read_all_quantified_peaks_file(self):
        try:
            all_quantified_peaks_path = self.search_task_path + 'AllQuantifiedPeaks.tsv'
            all_quantified_peaks_dataframe = pd.read_csv(all_quantified_peaks_path, sep='\t')
            return all_quantified_peaks_dataframe
        except Exception as e:
            print("Error reading AllQuantifiedPeaks.tsv file: ", e)

        return None
    
    def __read_results_file(self):
        try:
            results_path = self.search_task_path + 'results.txt'
            results_dataframe = pd.read_csv(results_path, sep='\t')
            return results_dataframe
        except Exception as e:
            print("Error reading results.txt file: ", e)

        return None 
       
    def __read_experimental_design_file(self):
        try:
            experimental_design_path = self.search_task_path + 'ExperimentalDesign.tsv'
            experimental_design_dataframe = pd.read_csv(experimental_design_path, sep='\t')
            return experimental_design_dataframe
        except Exception as e:
            print("Error reading ExperimentalDesign.tsv file: ", e)

        return None

class Plot:
    def __init__(self, bar_containers:List[matplotlib.container.BarContainer]) -> None:
        self.figure, self.ax = plt.subplots()
        self.bar_container_consumer(bar_containers)
        plt.legend()
        plt.plot()

    def bar_container_consumer(self, bar_containers:List[matplotlib.container.BarContainer]) -> matplotlib.figure.Figure:
        self.ax.add_container(bar_container for bar_container in bar_containers)

    @staticmethod
    def unique_full_sequence_bar_container(psms: List[PSM], label:str) -> matplotlib.container.BarContainer:
        number_of_distinct_full_sequences = {psm.full_sequence for psm in psms}
        return matplotlib.container.BarContainer(len(number_of_distinct_full_sequences), width=1, linewidth=0.7, label=label)

def charge_state_distribution(charge_states:List[int]) -> None:
    fig, ax = plt.subplots()

    ax.hist(charge_states)
    # ax.title = "Precursor Charge Distribution"
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize your results")

    parser.add_argument("search_task_directory",
                        metavar='search_task_dir',
                        type=str,
                        help="enter the search task directory path")

    args = parser.parse_args()
    directory = args.search_task_directory
    df = pd.read_csv(os.path.join(directory, "AllPeptides.psmtsv"), sep="\t")
    df = df.rename_axis("id").reset_index()
    df.to_sql(con=database, name=PSM.__tablename__, if_exists="append", index=True)
    
    results = session.query(PSM).all()
    charge_state_distribution([i.precursor_charge for i in results])