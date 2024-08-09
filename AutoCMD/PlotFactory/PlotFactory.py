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
from sqlalchemy import Column, Integer, String, Float, Boolean
import argparse

database = sa.create_engine('sqlite:///:memory:')
Session = sessionmaker(bind=database)
session = Session()
Base = declarative_base()

class PSM(Base):
    __tablename__ = "psms"

    id = Column(Integer, primary_key=True)
    file_name = Column("File Name", String(250))
    scan_number = Column("Scan Number", String(250))
    scan_retention_time = Column("Scan Retention Time", String(250))
    number_of_experimental_peaks = Column("Num Experimental Peaks", String(250))
    total_ion_current = Column("Total Ion Current", String(250))
    precursor_scan_number = Column("Precursor Scan Number", String(250))
    precursor_charge = Column("Precursor Charge", String(250))
    precursor_intensity = Column("Precursor Intensity", String(250))
    precursor_mz = Column("Precursor MZ", String(250))
    precursor_mass = Column("Precursor Mass", String(250))
    score = Column("Score", String(250))
    delta_score = Column("Delta Score", String(250))
    notch = Column("Notch", String(250))
    base_sequence = Column("Base Sequence", String(250))
    full_sequence = Column("Full Sequence", String(250))
    essential_sequence = Column("Essential Sequence", String(250))
    ambiguity_level = Column("Ambiguity Level", String(250))
    psm_count_unambigous_less_than_one_percent_q_value = Column("PSM Count (unambiguous, <0.01 q-value)", String(250))
    mods = Column("Mods", String(250))
    mods_chemical_formula = Column("Mods Chemical Formulas", String(250))
    mods_combined_chemical_formula = Column("Mods Combined Chemical Formula", String(250))
    number_of_variable_mods = Column("Num Variable Mods", String(250))
    missed_cleavages = Column("Missed Cleavages", String(250))
    peptide_monoisotopic_mass = Column("Peptide Monoisotopic Mass", String(250))
    mass_difference_in_daltons = Column("Mass Diff (Da)", String(250))
    mass_difference_in_ppm = Column("Mass Diff (ppm)", String(250))
    protein_accession = Column("Protein Accession", String(250))
    protein_name = Column("Protein Name", String(250))
    gene_name = Column("Gene Name", String(250))
    organism_name = Column("Organism Name", String(250))
    identified_sequence_variations = Column("Identified Sequence Variations", String(250))
    splice_sites = Column("Splice Sites", String(250))
    contaminant = Column("Contaminant", String(250))
    decoy = Column("Decoy", String(250))
    peptide_description = Column("Peptide Description", String(250))
    start_and_end_residues_in_protein = Column("Start and End Residues In Protein", String(250))
    previous_amino_acid = Column("Previous Amino Acid", String(250))
    next_amino_acid = Column("Next Amino Acid", String(250))
    theoretical_searched = Column("Theoreticals Searched", String(250))
    decoy_contaminant_target = Column("Decoy/Contaminant/Target", String(250))
    matched_ion_series = Column("Matched Ion Series", String(250))
    matched_ion_mass_to_charge_ratios = Column("Matched Ion Mass-To-Charge Ratios", String(250))
    matched_ion_mass_diff_in_daltons = Column("Matched Ion Mass Diff (Da)", String(250))
    matched_ion_mass_diff_in_ppm = Column("Matched Ion Mass Diff (Ppm)", String(250))
    matched_ion_intensities = Column("Matched Ion Intensities", String(250))
    matched_ion_counts = Column("Matched Ion Counts", String(250))
    normalized_spectral_angle = Column("Normalized Spectral Angle", String(250))
    localized_scores = Column("Localized Scores", String(250))
    improvement_possible = Column("Improvement Possible", String(250))
    cumulative_target = Column("Cumulative Target", String(250))
    cumulative_decoy = Column("Cumulative Decoy", String(250))
    q_value = Column("QValue", String(250))
    cumulative_target_notch = Column("Cumulative Target Notch", String(250))
    cumulative_decoy_notch = Column("Cumulative Decoy Notch", String(250))
    q_value_notch = Column("QValue Notch", String(250))
    posterior_error_probability = Column("PEP", String(250))
    posterior_error_probability_q_value = Column("PEP_QValue", String(250))


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

# class PSM:
#     def __init__(self, row: pd.Series) -> None:
#         self.file_name = row.get("File Name")
#         self.scan_number = row.get("Scan Number")
#         self.scan_retention_time = row.get("Scan Retention Time")
#         self.number_of_experimental_peaks = row.get("Num Experimental Peaks")
#         self.total_ion_current = row.get("Total Ion Current")
#         self.precursor_scan_number = row.get("Precursor Scan Number")
#         self.precursor_charge = row.get("Precursor Charge")
#         self.precursor_intensity = row.get("Precursor Intensity")
#         self.precursor_mz = row.get("Precursor MZ")
#         self.precursor_mass = row.get("Precursor Mass")
#         self.score = row.get("Score")
#         self.delta_score = row.get("Delta Score")
#         self.notch = row.get("Notch")
#         self.base_sequence = row.get("Base Sequence")
#         self.full_sequence = row.get("Full Sequence")
#         self.essential_sequence = row.get("Essential Sequence")
#         self.ambiguity_level = row.get("Ambiguity Level")
#         self.psm_count_unambigous_less_than_one_percent_q_value = row.get("PSM Count (unambiguous, <0.01 q-value)")
#         self.mods = row.get("Mods")
#         self.mods_chemical_formula = row.get("Mods Chemical Formula")
#         self.mods_combined_chemical_formula = row.get("Mods Combined Chemical Formula")
#         self.number_of_variable_mods = row.get("Num Variable Mods")
#         self.missed_cleavages = row.get("Missed Cleavages")
#         self.peptide_monoisotopic_mass = row.get("Peptide Monoisotopic Mass")
#         self.mass_difference_in_daltons = row.get("Mass Diff (Da)")
#         self.mass_difference_in_ppm = row.get("Mass Diff (ppm)")
#         self.protein_accession = row.get("Protein Accession")
#         self.protein_name = row.get("Protein Name")
#         self.gene_name = row.get("Gene Name")
#         self.organism_name = row.get("Organism Name")
#         self.identified_sequence_variations = row.get("Identified Sequence Variations")
#         self.splice_sites = row.get("Splice Sites")
#         self.contaminant = row.get("Contaminant")
#         self.decoy = row.get("Decoy")
#         self.peptide_description = row.get("Peptide Description")
#         self.start_and_end_residues_in_protein = row.get("Start and End Residues in Protein")
#         self.previous_amino_acid = row.get("Previous Amino Acid")
#         self.next_amino_acid = row.get("Next Amino Acid")
#         self.theoretical_searched = row.get("Theoretical Searched")
#         self.decoy_contaminant_target = row.get("Decoy/Contaminant/Target")
#         self.matched_ion_series = row.get("Matched Ion Series")
#         self.matched_ion_mass_to_charge_ratios = row.get("Matched Ion Mass-To-Charge Ratios")
#         self.matched_ion_mass_diff_in_daltons = row.get("Matched Ion Mass Diff (Da)")
#         self.matched_ion_mass_diff_in_ppm = row.get("Matched Ion Mass Diff (Ppm)")
#         self.matched_ion_intensities = row.get("Matched Ion Intensities")
#         self.matched_ion_counts = row.get("Matched Ion Counts")
#         self.normalized_spectral_angle = row.get("Normalized Spectral Angle")
#         self.localized_scores = row.get("Localized Scores")
#         self.improvement_possible = row.get("Improvement Possible")
#         self.cumulative_target = row.get("Cumulative Target")
#         self.cumulative_decoy = row.get("Cumulative Decoy")
#         self.q_value = row.get("QValue")
#         self.cumulative_target_notch = row.get("Cumulative Target Notch")
#         self.cumulative_decoy_notch = row.get("Cumulative Decoy Notch")
#         self.q_value_notch = row.get("QValue Notch")
#         self.posterior_error_probability = row.get("PEP")
#         self.posterior_error_probability_q_value = row.get("PEP_QValue")

#     @staticmethod
#     def get_all_psms(path:str) -> List:
#         all_psms_dataframe = pd.read_csv(os.path.join(path), sep="\t")
#         rows = all_psms_dataframe.to_dict(orient='records')
#         return [PSM(row) for row in rows]
    
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
    
    results = session.query(PSM).limit(10).all()
    print(results)