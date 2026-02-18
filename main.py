import os
import yaml
import pipeline_functions as steps
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

    # --- SETUP DIRECTORIES ---
    base_out = config['output_directory']
    
    # --- PIPELINE STEPS FOLDERS ---
    genome_cons_folder = os.path.join(base_out, "01_genomes", "conservation")
    genome_bio_folder = os.path.join(base_out, "01_genomes", "biosafety")
    blast_results_folder = os.path.join(base_out, "02_orthologs")
    msa_out_folder = os.path.join(base_out, "03_alignments")
    conserved_out_folder = os.path.join(base_out, "04_conserved_regions")
    biosafety_out_folder = os.path.join(base_out, "05_biosafety_check")
    final_report_folder = os.path.join(base_out, "06_final_report")

    # --- STEP 1: Download Genomes (Split Strategy) ---
    
    # 1.1 Download Conservation Species (Fungi)
    # Ensure Main Species is included in the conservation folder
    cons_species_list = config['conservation_targets'].copy()
    if config['main_species'] not in cons_species_list:
        cons_species_list.append(config['main_species'])
    
    print(f"\n[{steps.get_timestamp()}] === STEP 1a: Downloading Pathogens Genomes (Conservation) ===")
    steps.download_species_cds(
        {"species_list": cons_species_list}, 
        genome_cons_folder  # Salva na subpasta conservation
    )

    # 1.2 Download Biosafety Species (Plants/Insects/Human)
    print(f"\n[{steps.get_timestamp()}] === STEP 1b: Downloading Off-Target Genomes (Biosafety) ===")
    steps.download_species_cds(
        {"species_list": config['off_target_species']}, 
        genome_bio_folder   # Salva na subpasta biosafety
    )

    # --- STEP 2: Find Orthologs (BLAST) ---
    # Busca ortólogos apenas entre os fungos (pasta conservation)
    print(f"\n[{steps.get_timestamp()}] === STEP 2: Finding Orthologs ===")
    blast_config = {
        "main_species": config['main_species'],
        "target_species": config['conservation_targets'], 
        "evalue": config['blast_orthologs']['evalue'],
        "perc_identity": config['blast_orthologs']['perc_identity'],
        "min_coverage": config['blast_orthologs']['min_coverage'],
        "num_threads": config['blast_orthologs']['num_threads']
    }
    steps.find_orthologs_blast(blast_config, genome_cons_folder, blast_results_folder)

    # --- STEP 3: Multiple Sequence Alignment (MSA) ---
    # Alinha as sequências encontradas nos fungos
    print(f"\n[{steps.get_timestamp()}] === STEP 3: Multiple Sequence Alignment ===")
    msa_config = {
        "tool_msa_path": config['tools']['clustal_path'],
        "output_fmt": "fasta",
        "threads": config['blast_orthologs']['num_threads']
    }
    steps.perform_multiple_alignment(msa_config, blast_results_folder, genome_cons_folder, msa_out_folder)

    # --- STEP 4: Find Conserved Regions ---
    print(f"\n[{steps.get_timestamp()}] === STEP 4: Scanning for Conserved Regions ===")
    cons_config = {
        "min_conserved_length": config['conservation_scan']['min_conserved_length'],
        "conservation_threshold": config['conservation_scan']['conservation_threshold']
    }
    steps.find_conserved_regions(cons_config, msa_out_folder, conserved_out_folder)

    # --- STEP 5: Biosafety Check (Off-Target) ---
    # BLAST contra os organismos off-target (pasta biosafety)
    print(f"\n[{steps.get_timestamp()}] === STEP 5: Biosafety Check (Off-Target Analysis) ===")
    biosafety_config = {
        "species_list": config['off_target_species'],
        "risk_min_length": config['biosafety_risk']['risk_min_length'],
        "risk_min_identity": config['biosafety_risk']['risk_min_identity'],
        "evalue": config['biosafety_risk']['blast_evalue'],
        "word_size": config['biosafety_risk']['blast_word_size'],
        "num_threads": config['blast_orthologs']['num_threads']
    }
    steps.check_biosafety(biosafety_config, conserved_out_folder, genome_bio_folder, biosafety_out_folder)

    # --- STEP 6: Generate Annotation Report ---
    # Precisa do genoma principal para pegar a descrição (está na pasta conservation)
    print(f"\n[{steps.get_timestamp()}] === STEP 6: Annotating Candidates ===")
    report_config = {"main_species": config['main_species']}
    steps.generate_annotated_report(report_config, biosafety_out_folder, conserved_out_folder, genome_cons_folder, final_report_folder)

    # --- STEP 7: Apply Blacklist Filter ---
    print(f"\n[{steps.get_timestamp()}] === STEP 7: Functional Filtering (Blacklist) ===")
    annotated_csv = os.path.join(final_report_folder, "FINAL_CANDIDATES_ANNOTATED.csv")
    
    # Carrega a blacklist do config ou usa padrão
    blacklist = config.get('blacklist_terms', None)
    steps.apply_text_biosafety_filter(final_report_folder, final_report_folder, blacklist_list=blacklist)

    print(f"\n[{datetime.now()}] PIPELINE FINISHED SUCCESSFULLY.")

if __name__ == "__main__":
    main()
