#!/usr/bin/env python3
"""
This script consolidates multiple VCF files into one merged VCF file.
It replicates the bash script functionality:
  - Parses command-line options (verbose mode, output directory)
  - Prompts the user for the input VCF directory, reference genome build, and expected VCF version.
  - Offers metadata input (interactive mode or template generation) or skipping metadata.
  - Validates individual VCF files for header fileformat and reference genome.
  - Merges valid VCFs using vcfpy.
  - Appends metadata (if provided) to the merged VCF header.
  - Performs a final validation of the merged VCF.
  - Logs execution details to a log file.

Requirements:
  pip install vcfpy
"""

import os
import sys
import glob
import argparse
import logging
import datetime
import shutil
import re
from collections import OrderedDict

try:
    import vcfpy
except ImportError:
    sys.exit("Error: vcfpy package is required. Please install it with 'pip install vcfpy'.")

# Configure logging
LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(LOG_FILE)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)


def log_message(message, verbose=False):
    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message):
    log_message("CRITICAL ERROR: " + message)
    print("A critical error occurred. Check {} for details.".format(LOG_FILE))
    sys.exit(1)


def handle_non_critical_error(message):
    log_message("WARNING: " + message)
    print("Warning: " + message)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Consolidate multiple VCF files into a single merged VCF file."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose mode for detailed output."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Specify the output directory for the merged VCF file (defaults to the input directory if not provided).",
    )
    return parser.parse_args()

def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", "\"", "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace("\"", "\\\"")
    return f'"{escaped}"'


def build_sample_metadata_line(entries: "OrderedDict[str, str]") -> str:
    """Serialize an ordered mapping into a single ##SAMPLE metadata line."""

    id_value = entries.get("ID", "").strip()
    if not id_value:
        raise ValueError("Sample metadata must include a non-empty ID value.")

    parts = []
    for key, raw_value in entries.items():
        if raw_value is None:
            continue
        value = str(raw_value).strip()
        if not value and key != "ID":
            continue
        parts.append(f"{key}={_format_sample_metadata_value(value)}")

    serialized = ",".join(parts)
    return f"##SAMPLE=<{serialized}>"


def prompt_input(prompt_message, validation_regex=None, error_message="Invalid input."):
    while True:
        value = input(prompt_message + ": ").strip()
        if not value:
            print("Input cannot be empty. Please try again.")
            continue
        if validation_regex and not re.match(validation_regex, value):
            print(error_message)
            continue
        return value


def prompt_optional_input(prompt_message):
    return input(prompt_message + " (optional): ").strip()


def prompt_metadata_interactive():
    """
    Prompts the user to enter all metadata fields interactively.
    The metadata is written to final_metadata.txt as a single ##SAMPLE line.
    """
    print("Entering interactive metadata input mode...\n")
    
    sample_group_id = prompt_input(
        "Enter Sample Group Name/ID (used as the SampleGroup identifier)"
    )
    sample_description = prompt_optional_input("Enter Sample Group Description")

    # 1. REFERENCE_GENOME_BUILD
    ref_build = prompt_input("Enter Reference Genome Build (e.g., GRCh38)")
    
    # 2. GENOME_COMPLEXITY (Size, Ploidy, GC)
    genome_size = prompt_input("Enter Genome Size (e.g., 3.2Gb)", r"^[0-9]+(\.[0-9]+)?[Gg][Bb]$", "Genome Size must be in format like 3.2Gb")
    ploidy = prompt_input("Enter Genome Ploidy (e.g., 2)", r"^[1-9][0-9]*$", "Ploidy must be a positive integer")
    gc_content = prompt_input("Enter GC Content (e.g., 40%)", r"^[0-9]{1,3}%$", "GC Content must be a percentage (e.g., 40%)")
    
    # 3. SAMPLE_ORIGIN (Tissue, CollectionMethod, StorageConditions, TimeStored)
    tissue = prompt_input("Enter Tissue Type")
    collection_method = prompt_input("Enter Sample Collection Method (e.g., EDTA, RNAlater)")
    storage_cond = prompt_input("Enter Storage Conditions (e.g., -80C)")
    time_stored = prompt_input("Enter Time Stored (Days)", r"^[0-9]+$", "Time Stored must be a positive integer representing days")
    
    # 4. MATERIAL_TYPE (Type and IntegrityNumber)
    material_type = prompt_input("Enter Material Type (DNA, RNA, cDNA)")
    integrity_number = prompt_input("Enter Integrity Number (DIN/RIN)", r"^[0-9]+(\.[0-9]+)?$", "Integrity Number must be a numerical value")
    
    # 5. LIBRARY_CONSTRUCTION (Kit, Fragmentation, AdapterLigationEfficiency, PCRCycles)
    lib_kit = prompt_input("Enter Library Construction Kit")
    fragmentation = prompt_input("Enter Fragmentation Method (Mechanical, Enzymatic)")
    adapter_efficiency = prompt_input("Enter Adapter Ligation Efficiency (High, Medium, Low)")
    pcr_cycles = prompt_input("Enter PCR Cycles", r"^[0-9]+$", "PCR Cycles must be a positive integer")
    
    # 6. INPUT_QUALITY (A260_A280, A260_A230, DNAConcentration)
    a260_a280 = prompt_input("Enter Purity Ratio A260/A280", r"^[0-9]+(\.[0-9]+)?$", "Value must be numerical")
    a260_a230 = prompt_input("Enter Purity Ratio A260/A230", r"^[0-9]+(\.[0-9]+)?$", "Value must be numerical")
    dna_concentration = prompt_input("Enter DNA Concentration (e.g., ng/µL)", r"^[0-9]+(\.[0-9]+)?$", "Value must be numerical")
    
    # 7. ILLUMINA_SEQ (Instrument, FlowCell, ChannelMethod, ClusterDensity, QCSoftware)
    illumina_instrument = prompt_input("Enter Illumina Instrument")
    illumina_flowcell = prompt_input("Enter Illumina FlowCell")
    illumina_channel = prompt_input("Enter Illumina Channel Method")
    illumina_cluster = prompt_input("Enter Illumina Cluster Density")
    illumina_qc = prompt_input("Enter Illumina QC Software")
    
    # 8. ONT_SEQ (Instrument, FlowCellVersion, PoreType, BiasVoltage)
    ont_instrument = prompt_input("Enter ONT Instrument")
    ont_flowcell_version = prompt_input("Enter ONT FlowCell Version")
    ont_pore = prompt_input("Enter ONT Pore Type")
    ont_bias = prompt_input("Enter ONT Bias Voltage")
    
    # 9. PACBIO_SEQ (Instrument, SMRTCellType, ZMWDensity)
    pacbio_instrument = prompt_input("Enter PacBio Instrument")
    pacbio_smrtcell = prompt_input("Enter PacBio SMRT Cell Type")
    pacbio_zmwdensity = prompt_input("Enter PacBio ZMW Density")
    
    # 10. IONTORRENT_SEQ (Instrument, ChipType, pHCalibration, FlowOrder, IonSphereMetrics)
    iontorrent_instrument = prompt_input("Enter Ion Torrent Instrument")
    iontorrent_chip = prompt_input("Enter Ion Torrent Chip Type")
    iontorrent_ph = prompt_input("Enter Ion Torrent pH Calibration")
    iontorrent_floworder = prompt_input("Enter Ion Torrent Flow Order")
    iontorrent_ionsphere = prompt_input("Enter Ion Torrent IonSphere Metrics")
    
    # 11. PLATFORM_INDEPENDENT (Instrument, Pooling, SequencingKit, BaseCallingAlg, Q30, NormalizedCoverage, RunSpecificCalibration)
    platform_independent_instrument = prompt_input("Enter Platform-Independent Instrument")
    pooling = prompt_input("Enter Pooling Strategy")
    sequencing_kit = prompt_input("Enter Sequencing Kit")
    base_calling = prompt_input("Enter Base Calling Algorithm")
    q30 = prompt_input("Enter Q30")
    normalized_cov = prompt_input("Enter Normalized Coverage")
    run_calibration = prompt_input("Enter Run Specific Calibration")

    # 12. BIOINFO_ALIGNMENT (Tool, Software, Params, RefGenomeVers, RecalibrationSettings)
    align_tool = prompt_input("Enter Alignment Tool")
    align_software = prompt_input("Enter Alignment Software")
    align_params = prompt_input("Enter Alignment Parameters")
    align_refvers = prompt_input("Enter Reference Genome Version for Alignment")
    recalibration = prompt_input("Enter Recalibration Settings")

    # 13. BIOINFO_VARIANT_CALLING (Tool, Version, FilteringThresholds, DuplicateHandling, MQ)
    variant_tool = prompt_input("Enter Variant Calling Tool")
    variant_version = prompt_input("Enter Variant Calling Tool Version")
    filtering_thresholds = prompt_input("Enter Variant Filtering Thresholds")
    duplicate_handling = prompt_input("Enter Duplicate Handling Strategy")
    mq = prompt_input("Enter Mapping Quality (MQ)", r"^[0-9]+$", "Mapping Quality must be a positive integer")
    
    # 14. BIOINFO_POSTPROC (Normalization, Harmonization)
    normalization = prompt_input("Enter Post-Processing Normalization")
    harmonization = prompt_input("Enter Post-Processing Harmonization")
    
    # 15. SAMPLE_GROUP (LabName, LabMail, LabPhone, SampleCount, InclusionCriteria, ExclusionCriteria)
    lab_name = prompt_input("Enter Lab Name")
    lab_mail = prompt_input("Enter Lab Email")
    lab_phone = prompt_input("Enter Lab Phone")
    sample_count = prompt_input("Enter Sample Count", r"^[1-9][0-9]*$", "Sample Count must be a positive integer")
    inclusion = prompt_input("Enter Inclusion Criteria")
    exclusion = prompt_input("Enter Exclusion Criteria")
    
    # The key names follow the section_field convention expected by
    # MetaGap's parse_vcf_file importer so metadata values map onto the
    # appropriate SampleGroup-related models.
    metadata_entries = OrderedDict(
        [
            ("ID", sample_group_id),
            ("Description", sample_description),
            ("Reference_Genome_Build", ref_build),
            ("Reference_Genome_Build_Build_Name", ref_build),
            ("Genome_Complexity_Size", genome_size),
            ("Genome_Complexity_Ploidy", ploidy),
            ("Genome_Complexity_GC_Content", gc_content),
            ("Sample_Origin_Tissue", tissue),
            ("Sample_Origin_Collection_Method", collection_method),
            ("Sample_Origin_Storage_Conditions", storage_cond),
            ("Sample_Origin_Time_Stored", time_stored),
            ("Material_Type_Material_Type", material_type),
            ("Material_Type_Integrity_Number", integrity_number),
            ("Library_Construction_Kit", lib_kit),
            ("Library_Construction_Fragmentation", fragmentation),
            ("Library_Construction_Adapter_Ligation_Efficiency", adapter_efficiency),
            ("Library_Construction_PCR_Cycles", pcr_cycles),
            ("Input_Quality_A260_A280", a260_a280),
            ("Input_Quality_A260_A230", a260_a230),
            ("Input_Quality_DNA_Concentration", dna_concentration),
            ("Illumina_Seq_Instrument", illumina_instrument),
            ("Illumina_Seq_Flow_Cell", illumina_flowcell),
            ("Illumina_Seq_Channel_Method", illumina_channel),
            ("Illumina_Seq_Cluster_Density", illumina_cluster),
            ("Illumina_Seq_QC_Software", illumina_qc),
            ("ONT_Seq_Instrument", ont_instrument),
            ("ONT_Seq_Flow_Cell_Version", ont_flowcell_version),
            ("ONT_Seq_Pore_Type", ont_pore),
            ("ONT_Seq_Bias_Voltage", ont_bias),
            ("PacBio_Seq_Instrument", pacbio_instrument),
            ("PacBio_Seq_SMRT_Cell_Type", pacbio_smrtcell),
            ("PacBio_Seq_ZMW_Density", pacbio_zmwdensity),
            ("IonTorrent_Seq_Instrument", iontorrent_instrument),
            ("IonTorrent_Seq_Chip_Type", iontorrent_chip),
            ("IonTorrent_Seq_pH_Calibration", iontorrent_ph),
            ("IonTorrent_Seq_Flow_Order", iontorrent_floworder),
            ("IonTorrent_Seq_Ion_Sphere_Metrics", iontorrent_ionsphere),
            ("Platform_Independent", platform_independent_instrument),
            ("Platform_Independent_Pooling", pooling),
            ("Platform_Independent_Sequencing_Kit", sequencing_kit),
            ("Platform_Independent_Base_Calling_Alg", base_calling),
            ("Platform_Independent_Q30", q30),
            ("Platform_Independent_Normalized_Coverage", normalized_cov),
            ("Platform_Independent_Run_Specific_Calibration", run_calibration),
            ("Bioinfo_Alignment_Software", align_software),
            ("Bioinfo_Alignment_Params", align_params),
            ("Bioinfo_Alignment_Ref_Genome_Version", align_refvers),
            ("Bioinfo_Alignment_Recalibration_Settings", recalibration),
            ("Bioinfo_Variant_Calling_Tool", variant_tool),
            ("Bioinfo_Variant_Calling_Version", variant_version),
            ("Bioinfo_Variant_Calling_Filtering_Thresholds", filtering_thresholds),
            ("Bioinfo_Variant_Calling_Duplicate_Handling", duplicate_handling),
            ("Bioinfo_Variant_Calling_MQ", mq),
            ("Bioinfo_Post_Proc_Normalization", normalization),
            ("Bioinfo_Post_Proc_Harmonization", harmonization),
            ("Sample_Group_Lab_Name", lab_name),
            ("Sample_Group_Contact_Email", lab_mail),
            ("Sample_Group_Contact_Phone", lab_phone),
            ("Sample_Group_Total_Samples", sample_count),
            ("Sample_Group_Inclusion_Criteria", inclusion),
            ("Sample_Group_Exclusion_Criteria", exclusion),
        ]
    )
    metadata_line = build_sample_metadata_line(metadata_entries)

    with open("final_metadata.txt", "w") as mf:
        mf.write(metadata_line + "\n")

    print("Metadata has been successfully saved to final_metadata.txt")
    sys.exit(0)


def generate_template():
    # Template uses the same naming convention to remain import-compatible.
    template_entries = OrderedDict(
        [
            ("ID", "SampleGroupName"),
            ("Description", "Brief description of the cohort"),
            ("Reference_Genome_Build", "GRCh38"),
            ("Reference_Genome_Build_Build_Name", "GRCh38"),
            ("Genome_Complexity_Size", "3.2Gb"),
            ("Genome_Complexity_Ploidy", "2"),
            ("Genome_Complexity_GC_Content", "40%"),
            ("Sample_Origin_Tissue", "Blood"),
            ("Sample_Origin_Collection_Method", "EDTA"),
            ("Sample_Origin_Storage_Conditions", "-80C"),
            ("Sample_Origin_Time_Stored", "30"),
            ("Material_Type_Material_Type", "DNA"),
            ("Material_Type_Integrity_Number", "9.5"),
            ("Library_Construction_Kit", "Illumina TruSeq"),
            ("Library_Construction_Fragmentation", "Mechanical"),
            ("Library_Construction_Adapter_Ligation_Efficiency", "High"),
            ("Library_Construction_PCR_Cycles", "12"),
            ("Input_Quality_A260_A280", "1.8"),
            ("Input_Quality_A260_A230", "2.0"),
            ("Input_Quality_DNA_Concentration", "50"),
            ("Illumina_Seq_Instrument", "NovaSeq 6000"),
            ("Illumina_Seq_Flow_Cell", "S4"),
            ("Illumina_Seq_Channel_Method", "Dual"),
            ("Illumina_Seq_Cluster_Density", "250K/mm2"),
            ("Illumina_Seq_QC_Software", "BaseSpace"),
            ("ONT_Seq_Instrument", "MinION"),
            ("ONT_Seq_Flow_Cell_Version", "R9.4"),
            ("ONT_Seq_Pore_Type", "R9"),
            ("ONT_Seq_Bias_Voltage", "180mV"),
            ("PacBio_Seq_Instrument", "Sequel II"),
            ("PacBio_Seq_SMRT_Cell_Type", "8M"),
            ("PacBio_Seq_ZMW_Density", "50%"),
            ("IonTorrent_Seq_Instrument", "S5"),
            ("IonTorrent_Seq_Chip_Type", "530"),
            ("IonTorrent_Seq_pH_Calibration", "Standard"),
            ("IonTorrent_Seq_Flow_Order", "TACG"),
            ("IonTorrent_Seq_Ion_Sphere_Metrics", "Pass"),
            ("Platform_Independent", "NovaSeq 6000"),
            ("Platform_Independent_Pooling", "Multiplex"),
            ("Platform_Independent_Sequencing_Kit", "Universal"),
            ("Platform_Independent_Base_Calling_Alg", "Guppy"),
            ("Platform_Independent_Q30", "85%"),
            ("Platform_Independent_Normalized_Coverage", "30X"),
            ("Platform_Independent_Run_Specific_Calibration", "Complete"),
            ("Bioinfo_Alignment_Software", "BWA"),
            ("Bioinfo_Alignment_Params", "-t 8"),
            ("Bioinfo_Alignment_Ref_Genome_Version", "GRCh38"),
            ("Bioinfo_Alignment_Recalibration_Settings", "Default"),
            ("Bioinfo_Variant_Calling_Tool", "GATK"),
            ("Bioinfo_Variant_Calling_Version", "4.2"),
            ("Bioinfo_Variant_Calling_Filtering_Thresholds", "QD<2.0"),
            ("Bioinfo_Variant_Calling_Duplicate_Handling", "MarkDuplicates"),
            ("Bioinfo_Variant_Calling_MQ", "60"),
            ("Bioinfo_Post_Proc_Normalization", "LeftAlign"),
            ("Bioinfo_Post_Proc_Harmonization", "VEP"),
            ("Sample_Group_Lab_Name", "Genomics Lab"),
            ("Sample_Group_Contact_Email", "contact@genomicslab.org"),
            ("Sample_Group_Contact_Phone", "+1234567890"),
            ("Sample_Group_Total_Samples", "100"),
            ("Sample_Group_Inclusion_Criteria", "Age>18"),
            ("Sample_Group_Exclusion_Criteria", "None"),
        ]
    )

    template_line = build_sample_metadata_line(template_entries)

    with open("final_metadata_template.txt", "w") as tf:
        tf.write(template_line + "\n")
    print("Template final_metadata_template.txt created. Please fill it and rerun the script.")
    sys.exit(0)


def normalize_vcf_version(version_value):
    """Extract the numeric portion of the VCF version if present."""
    if not version_value:
        return version_value
    match = re.search(r"(\d+\.\d+)", version_value)
    if match:
        return match.group(1)
    return version_value.strip()


def find_first_vcf_with_header(input_dir, verbose=False):
    """Locate the first readable VCF file and return its header metadata."""
    candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    for file_path in candidates:
        try:
            reader = vcfpy.Reader.from_path(file_path)
        except Exception as e:
            handle_non_critical_error(
                f"Unable to open {file_path} for auto-detection: {str(e)}. Trying next file.",
            )
            continue

        header = reader.header
        fileformat = None
        reference = None
        for line in header.lines:
            if line.key == "fileformat" and fileformat is None:
                fileformat = line.value
            elif line.key == "reference" and reference is None:
                reference = line.value
            if fileformat and reference:
                break

        if fileformat and reference:
            log_message(
                f"Auto-detection using {file_path}: fileformat={fileformat}, reference={reference}",
                verbose,
            )
            return file_path, fileformat, reference

        handle_non_critical_error(
            f"Required metadata missing in {file_path} during auto-detection. Trying next file.",
        )

    return None, None, None


def validate_vcf(file_path, ref_genome, vcf_version, verbose=False):
    log_message(f"Validating file: {file_path}", verbose)
    if not os.path.isfile(file_path):
        handle_non_critical_error(f"File {file_path} does not exist. Skipping.")
        return False

    preprocessed_file = preprocess_vcf(file_path)
    try:
        reader = vcfpy.Reader.from_path(preprocessed_file)
    except Exception as e:
        handle_non_critical_error(f"Could not open {file_path}: {str(e)}. Skipping.")
        return False
    if preprocessed_file != file_path:
        os.remove(preprocessed_file)


    header = reader.header

    for line in header.lines:
        if isinstance(getattr(line, "key", None), str) and line.key.upper().startswith("GVCF"):
            handle_non_critical_error(
                f"{file_path} appears to be a gVCF (found header line {line.key}). Skipping."
            )
            return False

    fileformat = None
    for line in header.lines:
        if line.key == "fileformat":
            fileformat = line.value
            break

    expected_version = normalize_vcf_version(vcf_version)
    actual_version = normalize_vcf_version(fileformat)

    if fileformat is None or expected_version != actual_version:
        handle_non_critical_error(
            f"VCF version mismatch in {file_path}. Expected: {vcf_version}, Found: {fileformat}. Skipping.",
        )
        return False

    ref = None
    for line in header.lines:
        if line.key == "reference":
            ref = line.value
            break
    if ref is None or ref != ref_genome:
        handle_non_critical_error(f"Reference genome mismatch in {file_path}. Expected: {ref_genome}, Found: {ref}. Skipping.")
        return False


    try:
        from itertools import islice

        for record in islice(reader, 5):
            if any(
                (
                    (alt_value := getattr(alt, "value", "")) in {"<NON_REF>", "NON_REF"}
                    or str(alt) in {"<NON_REF>", "SymbolicAllele('NON_REF')"}
                )
                for alt in getattr(record, "ALT", [])
            ):
                handle_non_critical_error(
                    f"{file_path} appears to be a gVCF (found <NON_REF> ALT allele). Skipping."
                )
                return False

        _ = next(reader)
    except StopIteration:
        pass
    except Exception as e:
        handle_non_critical_error(f"Structural integrity check failed for {file_path}: {str(e)}. Skipping.")
        return False

    log_message(f"Validation passed for {file_path}", verbose)
    return True


def validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose=False):
    valid_vcfs = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    for file_path in glob.glob(os.path.join(input_dir, "*.vcf")):
        if validate_vcf(file_path, ref_genome, vcf_version, verbose):
            valid_vcfs.append(file_path)
        else:
            log_message(f"File {file_path} failed validation and is skipped.", verbose)
    if not valid_vcfs:
        handle_critical_error("No valid VCF files remain after validation. Aborting.")
    log_message("Validation completed. Valid VCF files: " + ", ".join(valid_vcfs), verbose)
    return valid_vcfs


def merge_vcfs(valid_files, output_dir, verbose=False):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    merged_filename = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    log_message("Merging VCF files...", verbose)
    try:
        reader0 = vcfpy.Reader.from_path(valid_files[0])
    except Exception as e:
        handle_critical_error("Failed to open the first valid VCF: " + str(e))
    header = reader0.header.copy()

    writer = vcfpy.Writer.from_path(merged_filename, header)
    for file_path in valid_files:
        log_message(f"Processing file: {file_path}", verbose)
        preprocessed_file = preprocess_vcf(file_path)
        try:
            reader = vcfpy.Reader.from_path(preprocessed_file)
            for record in reader:
                writer.write_record(record)
        except Exception as e:
            handle_non_critical_error(f"Error processing {file_path}: {str(e)}. Skipping remaining records.")
        if preprocessed_file != file_path:
            os.remove(preprocessed_file)

    writer.close()
    log_message(f"Merged VCF file created successfully: {merged_filename}", verbose)
    return merged_filename

def append_metadata_to_merged_vcf(merged_vcf, verbose=False):
    metadata_file = "final_metadata.txt"
    if not os.path.exists(metadata_file):
        log_message("No final_metadata.txt found. Skipping metadata append step.", verbose)
        return

    log_message("Appending final metadata to merged VCF header.", verbose)
    metadata_line = None
    with open(metadata_file, "r") as mf:
        for line in mf:
            stripped = line.strip()
            if stripped.startswith("##SAMPLE="):
                metadata_line = stripped
                break

    if metadata_line is None:
        log_message(
            "final_metadata.txt does not contain a ##SAMPLE line. Skipping metadata append.",
            verbose,
        )
        return

    if "<" not in metadata_line or not metadata_line.endswith(">"):
        log_message(
            "Malformed ##SAMPLE metadata detected; expected angle bracket encapsulation.",
            verbose,
        )
        return

    with open(merged_vcf, "r") as vf:
        merged_lines = vf.readlines()

    filtered_lines = [line for line in merged_lines if not line.startswith("##SAMPLE=")]

    insertion_index = next(
        (idx for idx, line in enumerate(filtered_lines) if line.startswith("#CHROM")),
        len(filtered_lines),
    )

    formatted_metadata_line = metadata_line
    if not formatted_metadata_line.endswith("\n"):
        formatted_metadata_line += "\n"

    filtered_lines.insert(insertion_index, formatted_metadata_line)

    with open(merged_vcf, "w") as vf:
        vf.writelines(filtered_lines)

    log_message("Metadata appended successfully.", verbose)

def preprocess_vcf(file_path):
    """
    Check if the VCF file uses spaces instead of tabs for the column header and data lines.
    If so, create a temporary file where the column header (#CHROM) and all subsequent lines 
    are converted to be tab-delimited. Lines starting with "##" (metadata) are left unchanged.
    Returns the path to the file to be used (original or temporary).
    """
    import re
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    modified = False
    new_lines = []
    header_found = False  # indicates when the column header has been encountered
    for line in lines:
        if line.startswith("##"):
            # Do not change metadata lines (they may contain spaces that are part of the value)
            new_lines.append(line)
        elif line.startswith("#"):
            # This is the column header line (e.g. "#CHROM ...")
            new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
        else:
            # Data lines: once header_found is True, convert spaces to tabs.
            if header_found:
                new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
                new_lines.append(new_line)
                if new_line != line:
                    modified = True
            else:
                new_lines.append(line)
    
    if modified:
        temp_file = file_path + ".tmp"
        with open(temp_file, 'w') as f:
            f.writelines(new_lines)
        return temp_file
    else:
        return file_path



def _info_field_requires_value(number):
    """Return True when the INFO definition expects an accompanying value."""

    if number is None:
        return False
    if isinstance(number, str):
        stripped = number.strip()
        if stripped == "":
            return False
        if stripped == "0":
            return False
        # symbolic numbers (., A, G, R) all expect a value when present
        if stripped in {".", "A", "G", "R"}:
            return True
        number = stripped
    try:
        return int(number) > 0
    except (TypeError, ValueError):
        # Unknown token, err on the side of requiring a value
        return True


def _has_null_value(value):
    """Return True if *value* or any nested value is ``None``."""

    if value is None:
        return True
    if isinstance(value, (list, tuple)):
        return any(_has_null_value(v) for v in value)
    return False


def validate_merged_vcf(merged_vcf, verbose=False):
    log_message(f"Starting validation of merged VCF: {merged_vcf}", verbose)
    if not os.path.isfile(merged_vcf):
        handle_critical_error(f"Merged VCF file {merged_vcf} does not exist.")
        
    try:
        reader = vcfpy.Reader.from_path(merged_vcf)
    except Exception as e:
        handle_non_critical_error(f"Could not open {merged_vcf}: {str(e)}. Skipping.")
        return False


    header = reader.header
    required_meta = ["fileformat", "reference"]
    for meta in required_meta:
        found = False
        for line in header.lines:
            if line.key == meta:
                found = True
                break
        if not found:
            handle_critical_error(f"Missing required meta-information: ##{meta} in {merged_vcf} header.")


    info_definitions = {}
    for line in header.lines:
        if isinstance(line, vcfpy.header.InfoHeaderLine):
            info_definitions[line.id] = line

    required_info_ids = {
        info_id for info_id, info_def in info_definitions.items() if _info_field_requires_value(info_def.number)
    }

    try:
        for record in reader:
            info_map = getattr(record, "INFO", None)
            record_label = f"{getattr(record, 'CHROM', '?')}:{getattr(record, 'POS', '?')}"

            if info_map is None:
                handle_critical_error(
                    f"Record {record_label} in {merged_vcf} is missing INFO data."
                )

            if not isinstance(info_map, dict):
                handle_critical_error(
                    f"Record {record_label} in {merged_vcf} has an unexpected INFO type: {type(info_map).__name__}."
                )

            if not info_map:
                missing_keys = sorted(required_info_ids)
            else:
                missing_keys = sorted(required_info_ids.difference(info_map.keys()))

            if missing_keys:
                handle_critical_error(
                    f"Record {record_label} in {merged_vcf} is missing required INFO fields: {', '.join(missing_keys)}."
                )

            null_keys = [key for key, value in info_map.items() if _has_null_value(value)]
            if null_keys:
                handle_critical_error(
                    f"Record {record_label} in {merged_vcf} has INFO fields with null values: {', '.join(sorted(null_keys))}."
                )

    except SystemExit:
        raise
    except Exception as exc:
        handle_critical_error(f"Error while parsing records in {merged_vcf}: {exc}")
    log_message(f"Validation completed successfully for merged VCF: {merged_vcf}", verbose)
    print(f"Validation completed successfully for merged VCF: {merged_vcf}")


def main():
    args = parse_arguments()
    verbose = args.verbose

    log_message("Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), verbose)

    input_dir = input("Enter the directory path containing VCF files: ").strip()
    while not os.path.isdir(input_dir):
        input_dir = input("Invalid directory. Please enter a valid directory path: ").strip()
    log_message("Input directory: " + input_dir, verbose)

    if args.output:
        output_dir = args.output
    else:
        temp_output_dir = input("Enter the output directory for the merged VCF file (Press Enter to use same as input directory): ").strip()
        output_dir = temp_output_dir if temp_output_dir else input_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    log_message("Output directory: " + output_dir, verbose)

    ref_genome = ""
    vcf_version = ""

    auto_choice = input(
        "Would you like to auto-detect the reference genome and VCF version from the first VCF header? (Y/N): "
    ).strip().upper()

    if auto_choice == "Y":
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(input_dir, verbose)
        if detected_file and detected_fileformat and detected_reference:
            detected_version = normalize_vcf_version(detected_fileformat)
            print("Auto-detected values:")
            print(f"  Source file: {detected_file}")
            print(f"  Reference genome: {detected_reference}")
            print(f"  VCF fileformat: {detected_fileformat} (version {detected_version})")
            confirm = input("Use these detected values? (Y/N): ").strip().upper()
            if confirm == "Y":
                ref_genome = detected_reference
                vcf_version = detected_version
            else:
                print("Proceeding with manual entry.")
        else:
            print("Auto-detection was unable to determine both reference and fileformat. Proceeding with manual entry.")

    while not ref_genome:
        user_input = input("Enter the reference genome build (e.g., GRCh38): ").strip()
        if user_input:
            ref_genome = user_input
        else:
            print("Reference genome build cannot be empty.")

    while not vcf_version:
        user_input = input("Enter the expected VCF version (e.g., 4.2): ").strip()
        if user_input:
            vcf_version = normalize_vcf_version(user_input)
        else:
            print("VCF version cannot be empty.")

    log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

    choice = input("Do you want to (I)nput metadata interactively, (T)emplate file, or (S)kip metadata? [I/T/S]: ").strip().upper()
    if choice == "I":
        prompt_metadata_interactive()
    elif choice == "T":
        if os.path.exists("final_metadata.txt"):
            print("Found final_metadata.txt. Using existing metadata file.")
        elif os.path.exists("final_metadata_template.txt"):
            use_template = input("final_metadata_template.txt exists. Do you want to use it as final_metadata.txt? (Y/N): ").strip().upper()
            if use_template == "Y":
                import shutil  # Ensure shutil is imported at the top
                shutil.copy("final_metadata_template.txt", "final_metadata.txt")
                print("Using final_metadata_template.txt as final_metadata.txt.")
            else:
                print("Please fill final_metadata_template.txt and rerun the script.")
                sys.exit(0)
        else:
            generate_template()
    elif choice == "S":
        print("Skipping metadata integration.")
    else:
        print("Invalid choice. No metadata will be integrated.")

    valid_files = validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose)
    merged_vcf = merge_vcfs(valid_files, output_dir, verbose)
    append_metadata_to_merged_vcf(merged_vcf, verbose)
    validate_merged_vcf(merged_vcf, verbose)

    print("----------------------------------------")
    print("Script Execution Summary:")
    print("Merged VCF File: " + merged_vcf)
    print("Log File: " + LOG_FILE)
    print("For detailed logs, refer to " + LOG_FILE)
    print("----------------------------------------")
    log_message("Script execution completed successfully.", verbose)


if __name__ == "__main__":
    main()
