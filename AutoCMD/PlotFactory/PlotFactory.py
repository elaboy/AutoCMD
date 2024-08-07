import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
import os

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
        self.file_name = row.get("File Name")
        self.scan_number = row.get("Scan Number")
        self.scan_retention_time = row.get("Scan Retention Time")
        self.number_of_experimental_peaks = row.get("Num Experimental Peaks")
        self.total_ion_current = row.get("Total Ion Current")
        self.precursor_scan_number = row.get("Precursor Scan Number")
        self.precursor_charge = row.get("Precursor Charge")
        self.precursor_intensity = row.get("Precursor Intensity")
        self.precursor_mz = row.get("Precursor MZ")
        self.precursor_mass = row.get("Precursor Mass")
        self.score = row.get("Score")
        self.delta_score = row.get("Delta Score")
        self.notch = row.get("Notch")
        self.base_sequence = row.get("Base Sequence")
        self.full_sequence = row.get("Full Sequence")
        self.essential_sequence = row.get("Essential Sequence")
        self.ambiguity_level = row.get("Ambiguity Level")
        self.psm_count_unambigous_less_than_one_percent_q_value = row.get("PSM Count (unambiguous, <0.01 q-value)")
        self.mods = row.get("Mods")
        self.mods_chemical_formula = row.get("Mods Chemical Formula")
        self.mods_combined_chemical_formula = row.get("Mods Combined Chemical Formula")
        self.number_of_variable_mods = row.get("Num Variable Mods")
        self.missed_cleavages = row.get("Missed Cleavages")
        self.peptide_monoisotopic_mass = row.get("Peptide Monoisotopic Mass")
        self.mass_difference_in_daltons = row.get("Mass Diff (Da)")
        self.mass_difference_in_ppm = row.get("Mass Diff (ppm)")
        self.protein_accession = row.get("Protein Accession")
        self.protein_name = row.get("Protein Name")
        self.gene_name = row.get("Gene Name")
        self.organism_name = row.get("Organism Name")
        self.identified_sequence_variations = row.get("Identified Sequence Variations")
        self.splice_sites = row.get("Splice Sites")
        self.contaminant = row.get("Contaminant")
        self.decoy = row.get("Decoy")
        self.peptide_description = row.get("Peptide Description")
        self.start_and_end_residues_in_protein = row.get("Start and End Residues in Protein")
        self.previous_amino_acid = row.get("Previous Amino Acid")
        self.next_amino_acid = row.get("Next Amino Acid")
        self.theoretical_searched = row.get("Theoretical Searched")
        self.decoy_contaminant_target = row.get("Decoy/Contaminant/Target")
        self.matched_ion_series = row.get("Matched Ion Series")
        self.matched_ion_mass_to_charge_ratios = row.get("Matched Ion Mass-To-Charge Ratios")
        self.matched_ion_mass_diff_in_daltons = row.get("Matched Ion Mass Diff (Da)")
        self.matched_ion_mass_diff_in_ppm = row.get("Matched Ion Mass Diff (Ppm)")
        self.matched_ion_intensities = row.get("Matched Ion Intensities")
        self.matched_ion_counts = row.get("Matched Ion Counts")
        self.normalized_spectral_angle = row.get("Normalized Spectral Angle")
        self.localized_scores = row.get("Localized Scores")
        self.improvement_possible = row.get("Improvement Possible")
        self.cumulative_target = row.get("Cumulative Target")
        self.cumulative_decoy = row.get("Cumulative Decoy")
        self.q_value = row.get("QValue")
        self.cumulative_target_notch = row.get("Cumulative Target Notch")
        self.cumulative_decoy_notch = row.get("Cumulative Decoy Notch")
        self.q_value_notch = row.get("QValue Notch")
        self.posterior_error_probability = row.get("PEP")
        self.posterior_error_probability_q_value = row.get("PEP_QValue")

    @staticmethod
    def get_all_psms(path:str):
        all_psms_dataframe = pd.read_csv(os.path.join(path), sep="\t")
        rows = all_psms_dataframe.to_dict(orient='records')
        return [PSM(row) for row in rows]
    
class QuantifiedPeak:
    def __init__(self, row:pd.Series):
        self.file_name = row.get("File Name")
        self.base_sequence = row.get("Base Sequence")
        self.full_sequence = row.get("Full Sequence")
        self.protein_group = row.get("Protein Group")
        self.peptide_monoisotopic_mass = row.get("Peptide Monoisotopic Mass")
        self.ms2_retention_time = row.get("MS2 Retention Time")
        self.precursor_charge = row.get("Precursor Charge")
        self.theoretical_mz = row.get("Theoretical MZ")
        self.peak_intensity = row.get("Peak Intensity")
        self.peak_retention_time_start = row.get("Peak RT Start")
        self.peak_retention_time_apex = row.get("Peak RT Apex")
        self.peak_retention_time_end = row.get("Peak RT End")
        self.peak_mz = row.get("Peak MZ")
        self.peak_charge = row.get("Peak Charge")
        self.number_charge_states_observed = row.get("Num Charge States Observed")
        self.peak_detection_type = row.get("Peak Detection Type")
        self.match_between_runs_score = row.get("MBR Score")
        self.psms_mapped = row.get("PSMs Mapped")
        self.base_sequences_mapped = row.get("Base Sequences Mapped")
        self.full_sequences_mapped = row.get("Full Sequences Mapped")
        self.peak_split_valley_retention_time = row.get("Peak Split Valley RT")
        self.peak_apex_mass_error_ppm = row.get("Peak Apex Mass Error (ppm)")
        self.match_between_runs_predicted_retention_time = row.get("MBR Predicted RT")

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
    