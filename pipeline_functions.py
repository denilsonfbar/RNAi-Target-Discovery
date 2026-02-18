import os
import sys
import shutil
import subprocess
import zipfile
from datetime import datetime
from glob import glob
import pandas as pd
from Bio import SeqIO
import json
from Bio import AlignIO


def find_key_recursive(obj, key_target):
    """
    Busca recursivamente por uma chave em um dicionário aninhado.
    """
    if isinstance(obj, dict):
        if key_target in obj:
            return obj[key_target]
        for key, value in obj.items():
            result = find_key_recursive(value, key_target)
            if result:
                return result
    elif isinstance(obj, list):
        for item in obj:
            result = find_key_recursive(item, key_target)
            if result:
                return result
    return None

def get_available_species_ncbi(genus):
    """
    Queries NCBI Datasets to list all species with available genomes for a given genus.
    """
    
    if shutil.which("datasets") is None:
        print(f"[{get_timestamp()}] Error: 'ncbi-datasets-cli' is not installed.")
        sys.exit(1)

    print(f"[{get_timestamp()}] Querying NCBI for available '{genus}' genomes...", end="", flush=True)

    cmd = [
        "datasets", "summary", "genome", "taxon", genus,
        "--as-json-lines",
        "--limit", "2000"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        unique_species = set()
        lines = result.stdout.splitlines()
        
        if not lines:
            print(" FAILED")
            print(f"[{get_timestamp()}] Error: NCBI returned no data.")
            return []

        # Debug var
        first_json_obj = None

        for line in lines:
            if not line.strip(): continue
            
            try:
                data = json.loads(line)
                
                if first_json_obj is None:
                    first_json_obj = data
                
                # --- CORREÇÃO AQUI ---
                # 1. Tenta a chave que vimos no seu Debug ('organism_name')
                full_name = find_key_recursive(data, "organism_name")
                
                # 2. Fallback para versões antigas ou diferentes ('sci_name')
                if not full_name:
                    full_name = find_key_recursive(data, "sci_name")
                
                # 3. Última tentativa ('scientific_name')
                if not full_name:
                    full_name = find_key_recursive(data, "scientific_name")
                
                if full_name:
                    parts = full_name.split()
                    # Garante binômio (Genero especie)
                    if len(parts) >= 2:
                        clean_name = f"{parts[0]} {parts[1]}"
                        unique_species.add(clean_name)
                        
            except json.JSONDecodeError:
                continue

        print(" OK")
        
        sorted_list = sorted(list(unique_species))
        
        if len(sorted_list) == 0:
            print(f"[{get_timestamp()}] Warning: 0 species found after parsing.")
            print(f"[{get_timestamp()}] DEBUG: Structure of the first JSON object keys:")
            if first_json_obj:
                print(list(first_json_obj.keys()))
                if 'organism' in first_json_obj:
                    print(f"Organism keys: {list(first_json_obj['organism'].keys())}")
            else:
                print("Could not parse JSON to inspect keys.")
        else:
            print(f"[{get_timestamp()}] Found {len(sorted_list)} unique species with genomes.")
        
        return sorted_list

    except subprocess.CalledProcessError as e:
        print(" FAILED")
        print(f"[{get_timestamp()}] NCBI Query Error: {e.stderr}")
        return []
    except Exception as e:
        print(f"\n[{get_timestamp()}] Unexpected Error: {e}")
        return []


def download_species_cds(config, output_folder):
    """
    Downloads the CDS (Coding DNA Sequences) for a list of species from NCBI.
    
    Logic:
    1. Iterates through the provided list of species.
    2. Tries to find a RefSeq assembly.
    3. If not found, checks for GenBank assembly.
    4. Prompts the user based on availability (RefSeq vs GenBank vs None).
    5. Downloads, extracts, and renames the file to {species_name}.fna.
    
    Args:
        config (dict): Configuration dictionary with "species_list" key.
        output_folder (str): Path to the destination folder.
    """
    
    species_list = config["species_list"]
    
    # 1. Setup Environment (Create folder once)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    # --- MAIN LOOP ---
    for species_name in species_list:
        
        # formatted species name for filenames (replace spaces with underscores)
        safe_name = species_name.replace(" ", "_")
        zip_path = os.path.join(output_folder, f"{safe_name}_temp.zip")
        temp_extract_dir = os.path.join(output_folder, f"temp_{safe_name}")

        # 2. Print Initial Log (Timestamp + Status)
        # end="" ensures the cursor stays on the same line to append "OK" later
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{current_time}] Downloading CDSs of {species_name}...", end="", flush=True)

        # 3. Check for RefSeq Availability using 'datasets summary'
        cmd_check_refseq = [
            "datasets", "summary", "genome", "taxon", species_name,
            "--assembly-source", "RefSeq",
            "--as-json-lines"
        ]
        
        # Run silently
        check_refseq = subprocess.run(cmd_check_refseq, capture_output=True, text=True)
        
        assembly_source = "RefSeq" # Default preference
        
        # If stdout is empty, RefSeq does not exist
        if not check_refseq.stdout.strip():
            
            # 3a. Check for GenBank Availability
            cmd_check_genbank = [
                "datasets", "summary", "genome", "taxon", species_name,
                "--as-json-lines"
            ]
            check_genbank = subprocess.run(cmd_check_genbank, capture_output=True, text=True)
            
            print("\n") # Move to next line to display the prompt clearly
            
            if check_genbank.stdout.strip():
                # Case: GenBank exists, but RefSeq does not
                print(f"   ⚠️  Warning: No RefSeq genome found for '{species_name}', but a GenBank assembly exists.")
                user_choice = input(f"   [1] Use GenBank assembly\n   [2] Skip this species\n   [3] Cancel pipeline\n   Select option: ")
                
                if user_choice == '1':
                    assembly_source = "GenBank"
                    # Reprint the status line so the final "OK" looks good
                    print(f"   Resuming download ({assembly_source})...", end="", flush=True)
                elif user_choice == '2':
                    print(f"   Skipping {species_name}.")
                    print("-" * 40)
                    continue # Skip to the next species in the loop
                else:
                    print("\nPipeline cancelled by user.")
                    sys.exit(1) # Stop the entire script
                    
            else:
                # Case: No genome found at all
                print(f"   ❌ Error: No genome found for '{species_name}' (neither RefSeq nor GenBank).")
                user_choice = input(f"   [1] Skip this species\n   [2] Cancel pipeline\n   Select option: ")
                
                if user_choice == '1':
                    print(f"   Skipping {species_name}.")
                    print("-" * 40)
                    continue # Skip to the next species
                else:
                    print("\nPipeline cancelled by user.")
                    sys.exit(1)

        # 4. Perform Download
        # Note: --limit flag removed as per previous debugging
        cmd_download = [
            "datasets", "download", "genome", "taxon", species_name,
            "--assembly-source", assembly_source,
            "--include", "cds",
            "--filename", zip_path
        ]
        
        try:
            # Run download (stderr piped to hide progress bar, ensuring clean output)
            subprocess.run(cmd_download, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            
            # 5. Process the Zip File
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(temp_extract_dir)

            # Find the CDS file inside the extracted structure
            # Pattern matches standard NCBI output: *_cds_from_genomic.fna
            found_files = glob(os.path.join(temp_extract_dir, "ncbi_dataset", "data", "*", "*cds_from_genomic.fna"))

            if not found_files:
                raise FileNotFoundError("CDS file not found inside the downloaded zip.")
                
            source_file = found_files[0]
            final_dest = os.path.join(output_folder, f"{safe_name}.fna")
            
            # Move and Rename
            shutil.move(source_file, final_dest)
            
            # 6. Success Output
            print(" OK")
            
        except subprocess.CalledProcessError:
            print(" FAILED")
            print(f"\n❌ Critical Error: Failed to download data for {species_name}.")
            # Option: sys.exit(1) to stop everything, or continue to try the next one?
            # Usually, if a download fails due to network, it's better to stop.
            sys.exit(1) 
        except Exception as e:
            print(" FAILED")
            print(f"\n❌ Error processing files: {e}")
            sys.exit(1)
        finally:
            # Cleanup temporary files
            if os.path.exists(zip_path):
                os.remove(zip_path)
            if os.path.exists(temp_extract_dir):
                shutil.rmtree(temp_extract_dir)


def get_timestamp():
    """Returns current date and time formatted as string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def count_genes_in_fasta(file_path):
    """Counts the number of sequences (headers starting with '>') in a FASTA file."""
    count = 0
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        print(f"[{get_timestamp()}] Error counting genes in {file_path}: {e}")
        return 0
    return count

def find_orthologs_blast(config, input_folder, output_folder):
    """
    Finds corresponding genes (orthologs) using BLASTn.
    
    Args:
        config (dict): Configuration dictionary with keys:
            - "main_species": Query species name.
            - "target_species_list": List of target species names.
            - "evalue": Maximum E-value to accept.
            - "perc_identity": Minimum percentage identity.
            - "min_coverage": Minimum query coverage percentage.
            - "num_threads": Number of CPU threads to use for BLAST.
        input_folder (str): Folder containing .fna files.
        output_folder (str): Folder to save resulting CSV files.
    """
    
    # 1. Setup
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    safe_main_name = config["main_species"].replace(" ", "_")
    query_file = os.path.join(input_folder, f"{safe_main_name}.fna")
    
    # Check Query File
    if not os.path.exists(query_file):
        print(f"[{get_timestamp()}] Error: Main species file not found: {query_file}")
        return {}

    # Count Query Genes
    total_query_genes = count_genes_in_fasta(query_file)
    
    print(f"[{get_timestamp()}] Starting BLAST Pipeline.")
    print(f"[{get_timestamp()}] Query Species: {config['main_species']} (Total Genes: {total_query_genes})")
    print(f"[{get_timestamp()}] Configuration: E-value < {config['evalue']}, Identity > {config['perc_identity']}%, Coverage > {config['min_coverage']}%")

    results_data = {}

    # 2. Iterate through Targets
    for target_species in config["target_species"]:
        if target_species == config["main_species"]:
            continue
            
        safe_target_name = target_species.replace(" ", "_")
        subject_file = os.path.join(input_folder, f"{safe_target_name}.fna")
        
        # Define DB name (stored in output folder to be reused)
        db_name = os.path.join(output_folder, f"blast_dbs/db_{safe_target_name}")
        
        # Define Output CSV path
        output_csv = os.path.join(output_folder, f"{safe_main_name}_vs_{safe_target_name}.csv")
        temp_blast_out = os.path.join(output_folder, f"temp_blast_{safe_target_name}.txt")

        # Check Target File
        if not os.path.exists(subject_file):
            print(f"[{get_timestamp()}] Skipping {target_species}: File not found.")
            continue

        # Count Target Genes
        total_target_genes = count_genes_in_fasta(subject_file)
        print(f"[{get_timestamp()}] Processing Target: {target_species} (Total Genes: {total_target_genes})")

        # 3. Check/Make BLAST DB
        if not os.path.exists(f"{db_name}.nhr"):
            print(f"[{get_timestamp()}] Building BLAST database for {target_species}...", end="", flush=True)
            cmd_makedb = [
                "makeblastdb", "-in", subject_file, "-dbtype", "nucl", 
                "-out", db_name, "-parse_seqids"
            ]
            subprocess.run(cmd_makedb, check=True, stdout=subprocess.DEVNULL)
            print(" OK")

        # 4. Run BLASTn
        # Columns: query_id, subject_id, identity, align_len, evalue, bitscore, query_len, subject_len
        cols = "6 qseqid sseqid pident length evalue bitscore qlen slen"
        
        print(f"[{get_timestamp()}] Running BLASTn...", end="", flush=True)
        
        cmd_blast = [
            "blastn",
            "-query", query_file,
            "-db", db_name,
            "-out", temp_blast_out,
            "-outfmt", cols,
            "-evalue", str(config['evalue']),
            "-perc_identity", str(config['perc_identity']),
            "-num_threads", str(config['num_threads']),
            "-max_target_seqs", "1" 
        ]
        
        try:
            subprocess.run(cmd_blast, check=True, stderr=subprocess.DEVNULL)
            print(" OK")
            
            # 5. Process Results
            if os.path.getsize(temp_blast_out) > 0:
                col_names = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qlen", "slen"]
                df = pd.read_csv(temp_blast_out, sep="\t", names=col_names)
          
                # Ordena pelo melhor score e remove duplicatas do mesmo gene Query
                # Isso garante que cada gene da espécie principal apareça apenas UMA vez na tabela
                df = df.sort_values("bitscore", ascending=False).drop_duplicates(subset=["qseqid"])      

                # Calculate Coverage
                df["coverage"] = (df["length"] / df["qlen"]) * 100
                
                # Filter
                df_filtered = df[df["coverage"] >= config['min_coverage']]
                
                # Save Final CSV
                df_filtered.to_csv(output_csv, index=False)
                
                # Store in dictionary
                results_data[target_species] = df_filtered
                
                # Log Stats
                print(f"[{get_timestamp()}] Raw Matches: {len(df)}")
                print(f"[{get_timestamp()}] Matches after Coverage Filter: {len(df_filtered)}")
                print(f"[{get_timestamp()}] Results saved to: {output_csv}")
                
            else:
                print(f"[{get_timestamp()}] No matches found (Empty BLAST output).")
                results_data[target_species] = pd.DataFrame()

            # Clean up temp file
            if os.path.exists(temp_blast_out):
                os.remove(temp_blast_out)

        except subprocess.CalledProcessError as e:
            print(f"[{get_timestamp()}] BLAST Execution Error: {e}")

    return results_data


def load_fasta_sequences(folder_path):
    """
    Loads ALL sequences from all .fna files in a folder into a dictionary.
    Structure: { 'Species_Name': { 'Gene_ID': SeqRecord } }
    """
    seq_db = {}
    files = [f for f in os.listdir(folder_path) if f.endswith(".fna")]
    
    print(f"[{get_timestamp()}] Loading sequences from {len(files)} genome files into memory...")
    
    for f in files:
        species_name = f.replace(".fna", "").replace("_", " ") # Revert filename to species name
        file_path = os.path.join(folder_path, f)
        
        # Dictionary for this species
        species_seqs = {}
        try:
            # We use SeqIO.to_dict for fast lookup
            # Note: We strip description to ensure ID matching matches BLAST output
            for record in SeqIO.parse(file_path, "fasta"):
                # Clean ID: sometimes BLAST returns IDs without 'lcl|' prefix
                clean_id = record.id.replace("lcl|", "")
                species_seqs[clean_id] = record
                # Also store original ID just in case
                species_seqs[record.id] = record
                
            seq_db[species_name] = species_seqs
            # Also store with underscore name for easier matching
            seq_db[f.replace(".fna", "")] = species_seqs
            
        except Exception as e:
            print(f"[{get_timestamp()}] Error reading {f}: {e}")
            
    print(f"[{get_timestamp()}] Sequences loaded. Memory ready.")
    return seq_db

def perform_multiple_alignment(config, blast_results_input_folder, genomes_input_folder, output_folder):
    """
    Parses BLAST CSVs, retrieves sequences, and runs Clustal Omega.
    """
    
    # 1. Setup
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Contar o total de espécies esperadas
    all_genome_files = [f for f in os.listdir(genomes_input_folder) if f.endswith(".fna")]
    total_expected_species = len(all_genome_files)
    
    print(f"[{get_timestamp()}] Analyzing {total_expected_species} species.")
    print(f"[{get_timestamp()}] Alignments will ONLY occur if a gene is present in ALL {total_expected_species} species.")

    # 2. Load Genomes (DNA)
    seq_database = load_fasta_sequences(genomes_input_folder)

    # 3. Group Hits by Query Gene
    gene_groups = {}
    
    csv_files = [f for f in os.listdir(blast_results_input_folder) if f.endswith(".csv")]
    print(f"[{get_timestamp()}] Parsing {len(csv_files)} BLAST result files...")

    for csv_file in csv_files:
        try:
            parts = csv_file.replace(".csv", "").split("_vs_")
            main_species_key = parts[0]
            target_species_key = parts[1]
        except IndexError:
            print(f"[{get_timestamp()}] Skipping invalid filename: {csv_file}")
            continue

        df = pd.read_csv(os.path.join(blast_results_input_folder, csv_file))
        
        for _, row in df.iterrows():
            q_id = row['qseqid']
            s_id = row['sseqid']
            
            if q_id not in gene_groups:
                gene_groups[q_id] = []
                gene_groups[q_id].append( {'species': main_species_key, 'gene': q_id} )
            
            gene_groups[q_id].append( {'species': target_species_key, 'gene': s_id} )

    # --- [ALTERAÇÃO 1] Filtragem Antes do Loop (Strict Mode) ---
    print(f"[{get_timestamp()}] Filtering candidates...")
    
    strict_groups = {}
    for q_id, members in gene_groups.items():
        # Verifica quantas espécies únicas estão no grupo
        unique_species = set(m['species'] for m in members)
        if len(unique_species) == total_expected_species:
            strict_groups[q_id] = members
            
    # Substitui o dicionário original pelo filtrado
    gene_groups = strict_groups
    total_groups = len(gene_groups)

    print(f"[{get_timestamp()}] Found {total_groups} VALID gene groups (present in all species).")
    
    if total_groups == 0:
        print(f"[{get_timestamp()}] ⚠️ No core genes found. Check BLAST thresholds.")
        return

    print(f"[{get_timestamp()}] Starting Clustal Omega alignment...")

    # 4. Run MSA for each group
    aligned_count = 0
    # skipped_count removido pois já filtramos tudo antes
    
    for i, (query_id, members) in enumerate(gene_groups.items()):
        
        # Progress Log
        print(f"   [{get_timestamp()}] Progress: {i}/{total_groups} aligned...", end="\r")

        # Collect SeqRecords
        sequences_to_align = []
        
        for member in members:
            sp = member['species']
            gene = member['gene']
            
            if sp in seq_database and gene in seq_database[sp]:
                record = seq_database[sp][gene]
                # Renomeia ID: >Citrus_clementina|XM_006423
                record.id = f"{sp}|{gene}"
                record.description = "" 
                sequences_to_align.append(record)

        # Se alguma sequência falhou em carregar (erro raro), pula
        if len(sequences_to_align) < total_expected_species:
            continue

        # --- [ALTERAÇÃO 2] Define extensão e nomes dos arquivos ---
        # Define extensão baseada na config (estética)
        ext = "fasta" if config['output_fmt'] in ['fasta', 'fa'] else "aln"
        
        # Usa ID único no temp para evitar conflitos
        temp_input = os.path.join(output_folder, f"temp_input_{query_id}.fasta")
        final_output = os.path.join(output_folder, f"msa_{query_id}.{ext}")
        
        # Write Temp Fasta
        SeqIO.write(sequences_to_align, temp_input, "fasta")
        
        # Build Clustal Command
        cmd_msa = [
            config['tool_msa_path'],
            "-i", temp_input,
            "-o", final_output,
            f"--outfmt={config['output_fmt']}",
            f"--threads={config['threads']}",
            "--force"
        ]
        
        try:
            subprocess.run(cmd_msa, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            aligned_count += 1
        except subprocess.CalledProcessError:
            print(f"\n[{get_timestamp()}] ❌ MSA Failed for group {query_id}")
        finally:
            if os.path.exists(temp_input):
                os.remove(temp_input)

    print(f"\n[{get_timestamp()}] Alignment Completed!")
    print(f"   Successfully Aligned: {aligned_count}")


def find_conserved_regions(config, msa_input_folder, output_folder):
    """
    Scans MSA files using a Sliding Window approach.
    FIX: Now includes Start-End positions in FASTA headers to avoid duplicates.
    """
    
    # 1. Setup
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    min_len = config.get("min_conserved_length", 60)
    threshold_pct = config.get("conservation_threshold", 100)
    threshold_ratio = threshold_pct / 100.0
    
    msa_files = glob(os.path.join(msa_input_folder, "*.fasta")) + glob(os.path.join(msa_input_folder, "*.aln"))
    
    print(f"[{get_timestamp()}] Scanning {len(msa_files)} alignments...")
    
    targets_found = []
    
    for msa_file in msa_files:
        try:
            filename = os.path.basename(msa_file)
            gene_id = filename.replace("msa_", "").split(".")[0]
            alignment = AlignIO.read(msa_file, "fasta")
            aln_len = alignment.get_alignment_length()
            
            # Step A: Consensus Mask
            consensus_mask = []
            consensus_seq = []
            for i in range(aln_len):
                col = alignment[:, i]
                unique_bases = set(col.upper())
                if '-' in unique_bases or 'N' in unique_bases:
                    consensus_mask.append(0)
                    consensus_seq.append('N')
                elif len(unique_bases) == 1:
                    consensus_mask.append(1)
                    consensus_seq.append(list(unique_bases)[0])
                else:
                    consensus_mask.append(0)
                    consensus_seq.append(list(unique_bases)[0])

            # Step B: Sliding Window
            i = 0
            while i <= aln_len - min_len:
                window = consensus_mask[i : i + min_len]
                matches = sum(window)
                identity = matches / min_len
                
                if identity >= threshold_ratio:
                    current_end = i + min_len 
                    while current_end < aln_len:
                        next_val = consensus_mask[current_end]
                        new_matches = matches + next_val
                        new_len = (current_end + 1) - i
                        new_identity = new_matches / new_len
                        
                        if new_identity >= threshold_ratio:
                            matches = new_matches
                            current_end += 1
                        else:
                            break
                    
                    block_seq = "".join(consensus_seq[i : current_end])
                    final_len = len(block_seq)
                    final_matches = sum(consensus_mask[i : current_end])
                    block_identity = (final_matches / final_len) * 100
                    
                    targets_found.append({
                        "gene_id": gene_id,
                        "start_pos": i,
                        "end_pos": current_end - 1,
                        "length": final_len,
                        "identity": round(block_identity, 2),
                        "sequence": block_seq,
                        "source_file": filename
                    })
                    i = current_end
                else:
                    i += 1

        except Exception as e:
            print(f"[{get_timestamp()}] Error processing {msa_file}: {e}")

    # 2. Save Results
    if targets_found:
        df = pd.DataFrame(targets_found)
        df = df.sort_values(["length", "identity"], ascending=[False, False])
        
        csv_path = os.path.join(output_folder, "conserved_targets.csv")
        df.to_csv(csv_path, index=False)
        
        fasta_path = os.path.join(output_folder, "conserved_targets.fasta")
        with open(fasta_path, "w") as f:
            for idx, row in df.iterrows():
                # --- CORREÇÃO AQUI: Adicionado '_posX-Y' para garantir unicidade ---
                header = f">{row['gene_id']}_len{row['length']}_id{row['identity']}_pos{row['start_pos']}-{row['end_pos']}"
                f.write(f"{header}\n{row['sequence']}\n")
        
        print(f"[{get_timestamp()}] Found {len(targets_found)} conserved regions.")
        print(f"[{get_timestamp()}] Results saved to: {csv_path}")
        return len(targets_found)
    else:
        print(f"[{get_timestamp()}] No regions found matching the criteria.")
        return 0


def check_biosafety(config, targets_input_folder, genomes_input_folder, output_folder):
    """
    Performs Off-Target analysis by BLASTing candidate targets against non-target genomes.
    
    Args:
        config (dict): Contains 'species_list', risk thresholds, and BLAST params.
        targets_input_folder (str): Path to the 'conserved_targets.fasta' file.
        genomes_input_folder (str): Folder where off-target genomes are stored.
        output_folder (str): Where to save biosafety reports.
    """
    
    # 1. Setup
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    print(f"[{get_timestamp()}] Starting Biosafety Check (Off-Target Analysis)...")
    
    # Extract BLAST parameters from config (with defaults)
    blast_evalue = str(config.get("evalue", "1e-5"))
    blast_word_size = str(config.get("word_size", "11"))
    blast_threads = str(config.get("num_threads", "4"))
    
    # Load all candidates first to track who survives
    candidates = []
    if os.path.exists(targets_input_folder):
        with open(os.path.join(targets_input_folder, "conserved_targets.fasta"), "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                candidates.append(record.id)
    else:
        print(f"[{get_timestamp()}] Error: Targets file not found: {targets_input_folder}")
        return []
            
    print(f"[{get_timestamp()}] Analyzing {len(candidates)} candidates against {len(config['species_list'])} genomes.")
    print(f"[{get_timestamp()}] BLAST Settings: E-value={blast_evalue}, Word Size={blast_word_size}, Threads={blast_threads}")

    # Dictionary to mark unsafe targets: { 'Gene_ID': ['Reason (Hit Citrus)', 'Reason (Hit Bee)'] }
    unsafe_targets = {}
    
    # 2. Iterate through Off-Target Species
    for species in config['species_list']:
        safe_sp_name = species.replace(" ", "_")
        genome_file = os.path.join(genomes_input_folder, f"{safe_sp_name}.fna")
        db_name = os.path.join(genomes_input_folder, f"db_{safe_sp_name}")
        blast_out = os.path.join(output_folder, f"blast_targets_vs_{safe_sp_name}.txt")
        
        # Check if genome exists
        if not os.path.exists(genome_file):
            print(f"[{get_timestamp()}] Skipping {species}: Genome file not found in {genomes_input_folder}.")
            continue
            
        # 3. Build DB (Reuse logic)
        if not os.path.exists(f"{db_name}.nhr"):
            print(f"[{get_timestamp()}] Building BLAST database for {species}...", end="", flush=True)
            cmd_makedb = ["makeblastdb", "-in", genome_file, "-dbtype", "nucl", "-out", db_name]
            subprocess.run(cmd_makedb, check=True, stdout=subprocess.DEVNULL)
            print(" OK")
            
        # 4. Run BLASTn
        cols = "6 qseqid sseqid pident length evalue"
        
        print(f"[{get_timestamp()}] Scanning against {species}...", end="", flush=True)
        
        cmd_blast = [
            "blastn",
            "-query", os.path.join(targets_input_folder, "conserved_targets.fasta"),
            "-db", db_name,
            "-out", blast_out,
            "-outfmt", cols,
            "-evalue", blast_evalue,
            "-word_size", blast_word_size,
            "-num_threads", blast_threads
        ]
        
        subprocess.run(cmd_blast, check=True)
        print(" OK")
        
        # 5. Analyze Hits (Filter "Risky" ones)
        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
            df = pd.read_csv(blast_out, sep="\t", names=["qseqid", "sseqid", "pident", "length", "evalue"])
            
            # CRITERIA FOR DANGER
            risk_len = config.get("risk_min_length", 19) 
            risk_id = config.get("risk_min_identity", 80.0)
            
            risky_hits = df[ (df['length'] >= risk_len) & (df['pident'] >= risk_id) ]
            
            unique_risky_ids = risky_hits['qseqid'].unique()
            
            if len(unique_risky_ids) > 0:
                print(f"      Found {len(unique_risky_ids)} potential off-target matches in {species}.")
                
                # Mark them as unsafe
                for gene_id in unique_risky_ids:
                    if gene_id not in unsafe_targets:
                        unsafe_targets[gene_id] = []
                    unsafe_targets[gene_id].append(species)
            else:
                 print(f"      No significant risks found against {species}.")
        else:
            print(f"      No hits found against {species}.")

    # 6. Generate Final Report
    safe_candidates = [gene for gene in candidates if gene not in unsafe_targets]
    
    # Save Safe List
    safe_csv = os.path.join(output_folder, "SAFE_targets_final.csv")
    pd.DataFrame(safe_candidates, columns=["gene_id"]).to_csv(safe_csv, index=False)
    
    # Save Unsafe List (with reasons)
    if unsafe_targets:
        unsafe_data = [{"gene_id": k, "risks": ", ".join(v)} for k, v in unsafe_targets.items()]
        unsafe_csv = os.path.join(output_folder, "UNSAFE_targets_rejected.csv")
        pd.DataFrame(unsafe_data).to_csv(unsafe_csv, index=False)
    
    print(f"\n[{get_timestamp()}] BIOSAFETY CHECK COMPLETE")
    print(f"   Total Candidates Checked: {len(candidates)}")
    print(f"   Rejected (Unsafe): {len(unsafe_targets)}")
    print(f"   APPROVED (Safe): {len(safe_candidates)}")
    print(f"   Final Approved List: {safe_csv}")
    
    return safe_candidates


def generate_annotated_report(config, final_targets_input_folder, conserved_targets_input_folder, genomes_input_folder, output_folder):
    """
    Merges safety data with functional annotations extracted from the original FASTA headers.
    FIX 1: Uses 'source_file' to reconstruct the FULL Gene ID for accurate description mapping.
    FIX 2: Uses 'composite_id' to accurately filter Safe targets from Conserved targets.
    """
    print(f"[{get_timestamp()}] 📝 Generating Final Annotated Report...")
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 1. Load Data
    try:
        # Load Safe IDs (IDs longos gerados na Etapa 6)
        safe_csv_path = os.path.join(final_targets_input_folder, "SAFE_targets_final.csv")
        safe_df = pd.read_csv(safe_csv_path)
        safe_ids_set = set(safe_df['gene_id'].astype(str).tolist())
        
        # Load Full Data (IDs curtos + Metadados da Etapa 5)
        conserved_csv_path = os.path.join(conserved_targets_input_folder, "conserved_targets.csv")
        conserved_df = pd.read_csv(conserved_csv_path)
        
        # --- FIX DE FILTRAGEM: Reconstruir o ID Composto ---
        # Recria o ID longo (ex: GeneA_len60_pos100-200) para cruzar com a lista SAFE
        conserved_df['composite_id'] = (
            conserved_df['gene_id'].astype(str) + 
            "_len" + conserved_df['length'].astype(str) + 
            "_id" + conserved_df['identity'].astype(str) + 
            "_pos" + conserved_df['start_pos'].astype(str) + "-" + conserved_df['end_pos'].astype(str)
        )
        
        # Filtra apenas as linhas cujo ID Composto está na lista de seguros
        final_df = conserved_df[conserved_df['composite_id'].isin(safe_ids_set)].copy()
        
        print(f"[{get_timestamp()}] Matched {len(final_df)} safe targets out of {len(conserved_df)} total regions.")
        
    except Exception as e:
        print(f"[{get_timestamp()}] Error loading CSV files: {e}")
        return

    # 2. Extract Annotations from Main Species FASTA
    main_species_file = os.path.join(genomes_input_folder, f"{config['main_species'].replace(' ', '_')}.fna")
    
    gene_descriptions = {}
    
    if not os.path.exists(main_species_file):
        print(f"[{get_timestamp()}] Warning: Main genome file not found at {main_species_file}. Annotations will be empty.")
    else:
        print(f"[{get_timestamp()}] Extracting gene descriptions from {os.path.basename(main_species_file)}...")
        
        try:
            with open(main_species_file, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    # O ID no arquivo FASTA é longo: lcl|NW_026850132.1_cds_XP_060400081.1_6396
                    full_id = record.id
                    
                    # Limpa a descrição (remove o ID repetido no início)
                    description = record.description.replace(full_id, "").strip()
                    
                    # Guarda no dicionário usando o ID COMPLETO como chave
                    gene_descriptions[full_id] = description
                    
                    # Opcional: Guarda também sem o prefixo lcl| para garantir
                    if "lcl|" in full_id:
                        gene_descriptions[full_id.replace("lcl|", "")] = description
                        
        except Exception as e:
             print(f"[{get_timestamp()}] Warning: Error parsing FASTA: {e}")

    # 3. Map Descriptions to DataFrame
    # --- FIX DE DESCRIÇÃO: Usar 'source_file' para obter o ID Completo ---
    
    def get_description_from_row(row):
        # O 'gene_id' na tabela está cortado (lcl|NW_...).
        # O 'source_file' tem o nome completo: msa_lcl|NW_026850132.1_cds_XP_060400081.1_6396.fasta
        
        filename = row['source_file']
        
        # 1. Remove extensão (.fasta, .aln) e prefixo (msa_)
        # os.path.splitext remove a última extensão
        base_name = os.path.splitext(filename)[0] 
        full_gene_id = base_name.replace("msa_", "")
        
        # 2. Tenta buscar este ID completo no dicionário
        desc = gene_descriptions.get(full_gene_id)
        
        if desc:
            return desc
        
        # 3. Tentativa de fallback (caso o ID tenha variações)
        return gene_descriptions.get(full_gene_id.replace("lcl|", ""), "No description found")

    # Aplica a função linha a linha
    final_df['description'] = final_df.apply(get_description_from_row, axis=1)

    # 4. Reorder Columns
    cols_order = ['gene_id', 'description', 'length', 'identity', 'start_pos', 'end_pos', 'sequence', 'source_file', 'composite_id']
    final_cols = [c for c in cols_order if c in final_df.columns]
    final_df = final_df[final_cols]

    # 5. Save Results
    csv_out = os.path.join(output_folder, "FINAL_CANDIDATES_ANNOTATED.csv")
    final_df.to_csv(csv_out, index=False)
    print(f"[{get_timestamp()}] 📊 Report saved: {csv_out}")
    
    print(f"[{get_timestamp()}] ✅ Done! You have {len(final_df)} annotated candidates ready for synthesis.")

def apply_text_biosafety_filter(input_folder, output_folder, blacklist_list=None):
    """
    Post-processing filter.
    Reads the annotated report and removes targets with descriptions matching specific
    housekeeping keywords (e.g., Ribosomes, Actin).
    """
    print(f"[{get_timestamp()}] 🧹 Starting Text-Based Functional Filter...")
    
    input_csv_path = os.path.join(input_folder, "FINAL_CANDIDATES_ANNOTATED.csv")
    if not os.path.exists(input_csv_path):
        print(f"[{get_timestamp()}] Error: Input file not found: {input_csv_path}")
        return

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    try:
        df = pd.read_csv(input_csv_path)
        print(f"[{get_timestamp()}] Loaded {len(df)} candidates for filtering.")
    except Exception as e:
        print(f"[{get_timestamp()}] Error reading CSV: {e}")
        return

    # Se o usuário não passou lista, usa uma padrão ou vazia
    if blacklist_list is None:
        risky_terms = ['ribosom', 'actin', 'tubulin'] # Padrão mínimo
    else:
        risky_terms = blacklist_list

    print(f"[{get_timestamp()}] Risky terms: {risky_terms}")
    
    def check_risk(desc):
        if pd.isna(desc): return None
        desc_lower = str(desc).lower()
        found_terms = []
        for term in risky_terms:
            if term in desc_lower:
                found_terms.append(term)
        return ", ".join(found_terms) if found_terms else None

    # Aplica o filtro
    df['risk_flag'] = df['description'].apply(check_risk)

    # Separa os DataFrames
    approved_df = df[df['risk_flag'].isnull()].copy()
    rejected_df = df[df['risk_flag'].notnull()].copy()
    
    # Remove a coluna auxiliar de flag do arquivo aprovado para ficar limpo
    if 'risk_flag' in approved_df.columns:
        approved_df = approved_df.drop(columns=['risk_flag'])

    # --- SAVE RESULTS ---
    
    # 1. TARGETS_GOLD_PREMIUM.csv (Os Aprovados)
    gold_csv = os.path.join(output_folder, "TARGETS_GOLD_PREMIUM.csv")
    approved_df.to_csv(gold_csv, index=False)
    
    # 2. TARGETS_REJECTED_FUNCTIONAL.csv (Os Rejeitados)
    rejected_csv = os.path.join(output_folder, "TARGETS_REJECTED_FUNCTIONAL.csv")
    rejected_df.to_csv(rejected_csv, index=False)

    print(f"[{get_timestamp()}] 🛡️  Filter Complete:")
    print(f"   - Input Candidates: {len(df)}")
    print(f"   - Rejected (Housekeeping): {len(rejected_df)}")
    print(f"   - ✅ APPROVED (GOLD): {len(approved_df)}")
    print(f"[{get_timestamp()}] 📂 Files saved to: {output_folder}")

    return rejected_df, approved_df
