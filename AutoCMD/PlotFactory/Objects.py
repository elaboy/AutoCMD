import sqlalchemy as sa
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
from sqlalchemy import Column, Integer, String, Float, Boolean
from typing import List, Self
from sqlalchemy.pool import StaticPool

engine = create_engine(
    "sqlite://", 
    connect_args={"check_same_thread": False}, 
    poolclass=StaticPool
)
Session = sessionmaker(bind=engine)
session = Session()

Base = declarative_base()
Base.metadata.create_all(engine)

class PSM(Base):
    __tablename__ = "psms"

    id = Column("id", Integer, primary_key=True)
    file_name = Column("File Name", String(250))
    scan_number = Column("Scan Number", Integer)
    scan_retention_time = Column("Scan Retention Time", Float)
    number_of_experimental_peaks = Column("Num Experimental Peaks", Integer)
    total_ion_current = Column("Total Ion Current", Float)
    precursor_scan_number = Column("Precursor Scan Number", Float)
    precursor_charge = Column("Precursor Charge", Integer)
    precursor_intensity = Column("Precursor Intensity", Float)
    precursor_mz = Column("Precursor MZ", Float)
    precursor_mass = Column("Precursor Mass", Float)
    score = Column("Score", Float)
    delta_score = Column("Delta Score", Float)
    notch = Column("Notch", Float)
    base_sequence = Column("Base Sequence", String(250))
    full_sequence = Column("Full Sequence", String(250))
    essential_sequence = Column("Essential Sequence", String(250))
    ambiguity_level = Column("Ambiguity Level", String(250))
    psm_count_unambigous_less_than_one_percent_q_value = Column("PSM Count (unambiguous, <0.01 q-value)", Float)
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
    matched_ion_counts = Column("Matched Ion Counts", Float)
    normalized_spectral_angle = Column("Normalized Spectral Angle", Float)
    localized_scores = Column("Localized Scores", Float)
    improvement_possible = Column("Improvement Possible", String(250))
    cumulative_target = Column("Cumulative Target", Float)
    cumulative_decoy = Column("Cumulative Decoy", Float)
    q_value = Column("QValue", Float)
    cumulative_target_notch = Column("Cumulative Target Notch", Float)
    cumulative_decoy_notch = Column("Cumulative Decoy Notch", Float)
    q_value_notch = Column("QValue Notch", Float)
    posterior_error_probability = Column("PEP", Float)
    posterior_error_probability_q_value = Column("PEP_QValue", Float)
    label = Column("Label", String(250))


class QuantifiedPeak(Base):
    __tablename__ = "quantified_peaks"

    id = Column(Integer, primary_key=True)
    file_name = Column("File Name", String(250))
    base_sequence = Column("Base Sequence", String(250))
    full_sequence = Column("Full Sequence", String(250))
    protein_group = Column("Protein Group", String(250))
    peptide_monoisotopic_mass = Column("Peptide Monoisotopic Mass", Float)
    ms2_retention_time = Column("MS2 Retention Time", Float)
    precursor_charge = Column("Precursor Charge", Integer)
    theoretical_mz = Column("Theoretical MZ", Float)
    peak_intensity = Column("Peak Intensity", Float)
    peak_retention_time_start = Column("Peak RT Start", Float)
    peak_retention_time_apex = Column("Peak RT Apex", Float)
    peak_retention_time_end = Column("Peak RT End", Float)
    peak_mz = Column("Peak MZ", Float)
    peak_charge = Column("Peak Charge", Integer)
    number_charge_states_observed = Column("Num Charge States Observed", Integer)
    peak_detection_type = Column("Peak Detection Type", String(250))
    match_between_runs_score = Column("MBR Score", Float)
    psms_mapped = Column("PSMs Mapped", Integer)
    base_sequences_mapped = Column("Base Sequences Mapped", Integer)
    full_sequences_mapped = Column("Full Sequences Mapped")
    peak_split_valley_retention_time = Column("Peak Split Valley RT", Float)
    peak_apex_mass_error_ppm = Column("Peak Apex Mass Error (ppm)", Float)
    match_between_runs_predicted_retention_time = Column("MBR Predicted RT", Float)
