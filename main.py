import os
import yaml
import pipeline_functions as steps  # Seu arquivo com as funções
from datetime import datetime

def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def main():
    # 1. Load Configuration
    config = load_config()
    
    print(f"--- RNAi Target Discovery Pipeline ---")
    print(f"Project: {config['project_name']}")
    print(f"Target: {config['main_species']}")
    print(f"Timestamp: {datetime.now()}")
    print("-" * 40)

    # Prepare Directories
    base_out = config['output_directory']
    genome_folder = os.path.join(base_out, "01_genomes")
    blast_db_folder = os.path.join(base_out, "02_blast_db")
    blast_results_folder = os.path.join(base_out, "03_orthologs")
    msa_out_folder = os.path.join(base_out, "04_alignments")
    conserved_out_folder = os.path.join(base_out, "05_conserved_regions")
    biosafety_out_folder = os.path.join(base_out, "06_biosafety_check")
    final_report_folder = os.path.join(base_out, "07_final_report")

    # --- STEP 1: Download Genomes ---
    # Combine lists for download (Fungi + Off-targets)
    all_species = list(set(config['conservation_targets'] + config['off_target_species']))
    
    # Create a temp config just for the download function
    download_config = {"species_list": all_species}
    steps.download_species_cds(download_config, genome_folder)

    # --- STEP 2: Find Orthologs (BLAST) ---
    blast_config = {
        "main_species": config['main_species'],
        "target_species": config['conservation_targets'], # Orthology check only among fungi
        "evalue": config['blast_orthologs']['evalue'],
        "perc_identity": config['blast_orthologs']['perc_identity'],
        "min_coverage": config['blast_orthologs']['min_coverage'],
        "num_threads": config['blast_orthologs']['num_threads']
    }
    steps.find_orthologs_blast(blast_config, genome_folder, blast_results_folder)

    # --- STEP 3: Multiple Sequence Alignment (MSA) ---
    msa_config = {
        "tool_msa_path": config['tools']['clustal_path'],
        "output_fmt": "fasta",
        "threads": config['blast_orthologs']['num_threads']
    }
    steps.perform_multiple_alignment(msa_config, blast_results_folder, genome_folder, msa_out_folder)

    # --- STEP 4: Find Conserved Regions ---
    cons_config = {
        "min_conserved_length": config['conservation_scan']['min_conserved_length'],
        "conservation_threshold": config['conservation_scan']['conservation_threshold']
    }
    steps.find_conserved_regions(cons_config, msa_out_folder, conserved_out_folder)

    # --- STEP 5: Biosafety Check (Off-Target) ---
    biosafety_config = {
        "species_list": config['off_target_species'],
        "risk_min_length": config['biosafety_risk']['risk_min_length'],
        "risk_min_identity": config['biosafety_risk']['risk_min_identity'],
        "evalue": config['biosafety_risk']['blast_evalue'],
        "word_size": config['biosafety_risk']['blast_word_size'],
        "num_threads": config['blast_orthologs']['num_threads']
    }
    steps.check_biosafety(biosafety_config, conserved_out_folder, genome_folder, biosafety_out_folder)

    # --- STEP 6: Generate Annotation Report ---
    report_config = {"main_species": config['main_species']}
    steps.generate_annotated_report(report_config, biosafety_out_folder, conserved_out_folder, genome_folder, final_report_folder)

    # --- STEP 7: Apply Blacklist Filter ---
    # Need to pass the blacklist from config to the function (Requires small update in steps.py or pass as arg)
    # For now, let's assume the function reads the hardcoded list or we update it.
    # BEST PRACTICE: Update apply_text_biosafety_filter to accept a list.
    
    annotated_csv = os.path.join(final_report_folder, "FINAL_CANDIDATES_ANNOTATED.csv")
    steps.apply_text_biosafety_filter(annotated_csv, final_report_folder, config['blacklist_terms'])

    print(f"\n[{datetime.now()}] PIPELINE FINISHED SUCCESSFULLY.")

if __name__ == "__main__":
    main()
