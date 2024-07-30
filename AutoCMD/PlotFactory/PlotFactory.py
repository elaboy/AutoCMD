import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Results:
    def __init__(self, search_task_path:str) -> None:
        self.search_task_path = search_task_path
        self.all_psms_dataframe = self.__read_all_psms_file()
        self.all_peptides_dataframe = self.__read_all_peptides_file()
        self.all_quantified_peaks_dataframe = self.__read_all_quantified_peaks_file()
        self.results_dataframe = self.__read_results_file()
        self.experimental_design_dataframe = self.__read_experimental_design_file()

    def __read_all_psms_file(self):
        all_psms_path = self.search_task_path + 'AllPSMs.psmtsv'
        all_psms_dataframe = pd.read_csv(all_psms_path, sep='\t')
        return all_psms_dataframe
    
    def __read_all_peptides_file(self):
        all_peptides_path = self.search_task_path + 'AllPeptides.psmtsv'
        all_peptides_dataframe = pd.read_csv(all_peptides_path, sep='\t')
        return all_peptides_dataframe
    
    def __read_all_quantified_peaks_file(self):
        all_quantified_peaks_path = self.search_task_path + 'AllQuantifiedPeaks.tsv'
        all_quantified_peaks_dataframe = pd.read_csv(all_quantified_peaks_path, sep='\t')
        return all_quantified_peaks_dataframe
    
    def __read_results_file(self):
        results_path = self.search_task_path + 'results.txt'
        results_dataframe = pd.read_csv(results_path, sep='\t')
        return results_dataframe

    def __read_experimental_design_file(self):
        experimental_design_path = self.search_task_path + 'ExperimentalDesign.tsv'
        experimental_design_dataframe = pd.read_csv(experimental_design_path, sep='\t')
        return experimental_design_dataframe

# this class will hold all the plotting functions, where
# each function will make one figure, and a following function will join them
class Plots:
