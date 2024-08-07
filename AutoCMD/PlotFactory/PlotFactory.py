import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum

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

class PSM:
    def __init__(self, row: pd.Series) -> None:
        self.file_name = row["File Name"]
        self.scan_number = row["Scan Number"]
        self.scan_retention_time = row["Scan Retention Time"]
        self.number_of_experimental_peaks = row["Num Experimental Peaks"]
        self.total_ion_current = row["Total Ion Current"]
        self.precursor_scan_number = row["Precursor Scan Number"]
        self.precursor_charge = row["Precursor Charge"]
        self.precursor_intensity = row["Precursor Intensity"]
        self.precursor_mz = row["Precursor MZ"]
        self.precursor_mass = row["Precursor Mass"]
        self.score = row["Score"]
        self.delta_score = row["Delta Score"]
        self.notch = row["Notch"]
        self.base_sequence = row["Base Sequence"]
        self.full_sequence = row["Full Sequence"]
        self.essential_sequence = row["Essential Sequence"]
        self.ambiguity_level = row["Ambiguity Level"]
        self.psm_count_unambigous_less_than_one_percent_q_value = row["PSM Count (unambiguous, <0.01 q-value)"]
        self.mods = row["Mods"]
        self.mods_chemical_formula = row["Mods Chemical Formula"]
        self.mods_combined_chemical_formula = row["Mods Combined Chemical Formula"]
        self.number_of_variable_mods = row["Num Variable Mods"]
        self.missed_cleavages = row["Missed Cleavages"]
        self.peptide_monoisotopic_mass = row["Peptide Monoisotopic Mass"]
        self.mass_difference_in_daltons = row["Mass Diff (Da)"]
        self.mass_difference_in_ppm = row["Mass Diff (ppm)"]
        self.protein_accession = row["Protein Accession"]
        self.protein_name = row["Protein Name"]
        self.gene_name = row["Gene Name"]
        self.organism_name = row["Organism Name"]
        self.identified_sequence_variations = row["Identified Sequence Variations"]
        self.splice_sites = row["Splice Sites"]
        self.contaminant = row["Contaminant"]
        self.decoy = row["Decoy"]
        self.peptide_description = row["Peptide Description"]
        self.start_and_end_residues_in_protein = row["Start and End Residues in Protein"]
        self.previous_amino_acid = row["Previous Amino Acid"]
        self.next_amino_acid = row["Next Amino Acid"]
        self.theoretical_searched = row["Theoretical Searched"]
        self.decoy_contaminant_target = row["Decoy/Contaminant/Target"]
        self.matched_ion_series = row["Matched Ion Series"]
        self.matched_ion_mass_to_charge_ratios = row["Matched Ion Mass-To-Charge Ratios"]
        self.matched_ion_mass_diff_in_daltons = row["Matched Ion Mass Diff (Da)"]
        self.matched_ion_mass_diff_in_ppm = row["Matched Ion Mass Diff (Ppm)"]
        self.matched_ion_intensities = row["Matched Ion Intensities"]
        self.matched_ion_counts = row["Matched Ion Counts"]
        self.normalized_spectral_angle = row["Normalized Spectral Angle"]
        self.localized_scores = row["Localized Scores"]
        self.improvement_possible = row["Improvement Possible"]
        self.cumulative_target = row["Cumulative Target"]
        self.cumulative_decoy = row["Cumulative Decoy"]
        self.q_value = row["QValue"]
        self.cumulative_target_notch = row["Cumulative Target Notch"]
        self.cumulative_decoy_notch = row["Cumulative Decoy Notch"]
        self.q_value_notch = row["QValue Notch"]
        self.posterior_error_probability = row["PEP"]
        self.posterior_error_probability_q_value = row["PEP_QValue"]

class QuantifiedPeak:
    def __init__(self, row:pd.Series):
        self.file_name = row["File Name"]
        self.base_sequence = row["Base Sequence"]
        self.full_sequence = row["Full Sequence"]
        self.protein_group = row["Protein Group"]
        self.peptide_monoisotopic_mass = row["Peptide Monoisotopic Mass"]
        self.ms2_retention_time = row["MS2 Retention Time"]
        self.precursor_charge = row["Precursor Charge"]
        self.theoretical_mz = row["Theoretical MZ"]
        self.peak_intensity = row["Peak Intensity"]
        self.peak_retention_time_start = row["Peak RT Start"]
        self.peak_retention_time_apex = row["Peak RT Apex"]
        self.peak_retention_time_end = row["Peak RT End"]
        self.peak_mz = row["Peak MZ"]
        self.peak_charge = row["Peak Charge"]
        self.number_charge_states_observed = row["Num Charge States Observed"]
        self.peak_detection_type = row["Peak Detection Type"]
        self.match_between_runs_score = row["MBR Score"]
        self.psms_mapped = row["PSMs Mapped"]
        self.base_sequences_mapped = row["Base Sequences Mapped"]
        self.full_sequences_mapped = row["Full Sequences Mapped"]
        self.peak_split_valley_retention_time = row["Peak Split Valley RT"]
        self.peak_apex_mass_error_ppm = row["Peak Apex Mass Error (ppm)"]
        self.match_between_runs_predicted_retention_time = row["MBR Predicted RT"]

class QuantifiedPeptide:
    def __init__(self, row:pd.Series):
        self.sequence = row["Sequence"]
        self.base_sequence = row["Base Sequence"]
        self.protein_groups = row["Protein Groups"]
        self.gene_names = row["Gene Names"]
        self.organism = row["Organism"]
        self.file_intensity = row[5]
        self.detection_type = row[6]

class QuantifiedProteinGroup:
    def __init(self, row:pd.Series):
        self.protein_accession = row["Protein Accession"]
        self.gene = row["Gene"]
        self.organism = row["Organism"]
        self.protein_full_name = row["Protein Full Name"]
        self.protein_unmodified_mass = row["Protein Unmodified Mass"]
        self.number_of_peptides_in_group = row["Number of Peptides in Group"]
        self.unique_peptides = row["Unique Peptides"]
        self.shared_peptides = row["Shared Peptides"]
        self.number_of_unique_peptides = row["Number of Unique Peptides"]
        self.sequence_coverage_fraction = row["Sequence Coverage Fraction"]
        self.sequence_coverage = row["Sequence Coverage"]
        self.sequence_coverage_with_mods = row["Sequence Coverage with Mods"]
        self.fragment_sequence_coverage = row["Fragment Sequence Coverage"]
        self.modification_info_list = row["Modification Info List"]
        self.file_intensity = row[14]
        self.number_of_psms = row["Number of PSMs"]
        self.protein_decoy_contaminant_target = row["Protein Decoy/Contaminant/Target"]
        self.protein_cumulative_target = row["Protein Cumulative Target"]
        self.protein_cumulative_decoy = row["Protein Cumulative Decoy"]
        self.protein_q_value = row["Protein QValue"]
        self.best_peptide_score = row["Best Peptide Score"]
        self.best_peptide_notch_q_value = row["Best Peptide Notch QValue"]

        
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