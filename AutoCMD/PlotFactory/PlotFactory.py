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
    @staticmethod
    def unique_base_sequences(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame) -> plt.Figure:
        number_of_distinct_base_sequences_df1 = dataframe1["Base Sequence"].nunique()
        number_of_distinct_base_sequences_df2 = dataframe2["Base Sequence"].nunique()

        fig, ax = plt.subplots()

        ax.bar(1, number_of_distinct_base_sequences_df1, width = 1, linewidth = 0.7, label="SSRCalc3")
        ax.bar(4, number_of_distinct_base_sequences_df2, width = 1, linewidth = 0.7, label="Chronologer")

        ax.set(xlim=(0, 8), xticks=np.arange(1, 6))

        #removes the X Axis ticks
        plt.xticks([])

        ax.set_title("Unique Base Sequences in All Peptides File")

        #dataframe1 annotation
        ax.annotate("n = "+ ("{:,}").format(number_of_distinct_base_sequences_df1),
                    xy=(1.3, number_of_distinct_base_sequences_df1-100000),
                    xytext=(2, number_of_distinct_base_sequences_df1),
                    arrowprops=dict(facecolor='red', width=2, headwidth=12))
        
        #dataframe two annotation
        ax.annotate("n = "+ ("{:,}").format(number_of_distinct_base_sequences_df2),
                    xy=(4.3, number_of_distinct_base_sequences_df2-100000),
                    xytext=(6, number_of_distinct_base_sequences_df2),
                    arrowprops=dict(facecolor='red', width=2, headwidth=12))

        ax.legend()

        return fig
    
    @staticmethod
    def unique_full_sequences(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame) -> plt.Figure:
        number_of_distinct_base_sequences_df1 = dataframe1["Full Sequence"].nunique()
        number_of_distinct_base_sequences_df2 = dataframe2["Full Sequence"].nunique()

        fig, ax = plt.subplots()

        ax.bar(1, number_of_distinct_base_sequences_df1, width = 1, linewidth = 0.7, label="SSRCalc3")
        ax.bar(4, number_of_distinct_base_sequences_df2, width = 1, linewidth = 0.7, label="Chronologer")

        ax.set(xlim=(0, 8), xticks=np.arange(1, 6))

        #removes the X Axis ticks
        plt.xticks([])

        ax.set_title("Unique Full Sequences in All Peptides File")

        #dataframe1 annotation
        ax.annotate("n = "+ ("{:,}").format(number_of_distinct_base_sequences_df1),
                    xy=(1.3, number_of_distinct_base_sequences_df1-100000),
                    xytext=(2, number_of_distinct_base_sequences_df1),
                    arrowprops=dict(facecolor='red', width=2, headwidth=12))
        
        #dataframe two annotation
        ax.annotate("n = "+ ("{:,}").format(number_of_distinct_base_sequences_df2),
                    xy=(4.3, number_of_distinct_base_sequences_df2-100000),
                    xytext=(6, number_of_distinct_base_sequences_df2),
                    arrowprops=dict(facecolor='red', width=2, headwidth=12))

        ax.legend()

        return fig