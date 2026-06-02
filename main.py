import pandas as pd
import pickle
import subprocess
import sys
import os
import argparse
import shutil
from pathlib import Path

def setup_abricate_databases():
    """Setup abricate databases if not already done"""
    # Check if VF_custom database is already set up by looking for BLAST index files
    db_path = "./abricate_reproducible/db/VF_custom/sequences"
    index_files = [f"{db_path}.nhr", f"{db_path}.nin", f"{db_path}.nsq"]
    
    # If all index files exist, databases are already set up
    if all(os.path.exists(f) for f in index_files):
        return True
    
    # Otherwise, run setup
    cmd = [
        "./abricate_reproducible/bin/abricate",
        "--setupdb"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print("Abricate databases setup completed successfully", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error setting up abricate databases: {e}", file=sys.stderr)
        return False

def run_blastn(input_fasta):
    """Run blastn command and return the output file path"""
    db_path = "ecoli/e.coli_pathotype.fa"
    output_file = f"{input_fasta}.blastn"
    
    cmd = [
        "blastn",
        "-db", db_path,
        "-query", input_fasta,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen",
        "-perc_identity", "90"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error running blastn: {e}", file=sys.stderr)
        return None

def run_abricate(input_fasta):
    """Run abricate commands"""
    # Create temp directory if it doesn't exist
    os.makedirs("temp", exist_ok=True)
    
    # First abricate command
    cmd1 = [
        "./abricate_reproducible/bin/abricate",
        "--nopath", "--quiet", "--db", "VF_custom",
        "--minid", "90", "--mincov", "60",
        input_fasta
    ]
    
    try:
        with open("temp/abr.abricate", "w") as f:
            subprocess.run(cmd1, stdout=f, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running first abricate command: {e}", file=sys.stderr)
        return False
    
    # Second command - abricate parser
    cmd2 = [
        "python3", "./abricate_reproducible/bin/abricate_parser.py",
        "temp/abr.abricate", "temp/abr.summarized.abricate", "\t"
    ]
    
    try:
        subprocess.run(cmd2, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running abricate parser: {e}", file=sys.stderr)
        return False
    
    # Third command - abricate summary
    cmd3 = [
        "./abricate_reproducible/bin/abricate",
        "--summary", "temp/abr.summarized.abricate"
    ]
    
    try:
        with open("temp/abr.tsv", "w") as f:
            subprocess.run(cmd3, stdout=f, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running abricate summary: {e}", file=sys.stderr)
        return False

def cleanup_files(blast_file):
    """Clean up temporary files"""
    try:
        # Remove blastn output file
        if os.path.exists(blast_file):
            os.remove(blast_file)
        
        # Remove only specific temporary files, not the entire temp directory
        temp_files = [
            "temp/abr.abricate",
            "temp/abr.summarized.abricate", 
            "temp/abr.tsv"
        ]
        
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                
    except Exception as e:
        # Don't fail if cleanup fails
        print(f"Warning: Could not clean up temporary files: {e}", file=sys.stderr)

pd.set_option('future.no_silent_downcasting', True)

def main():
    parser = argparse.ArgumentParser(description="RF prediction pipeline for pathogenicity assessment")
    parser.add_argument("input_fasta", help="Input FASTA file path")
    args = parser.parse_args()
    
    # Hardcoded paths
    model_path = "final_RF.pkl"
    columns_path = "test_df.csv"
    
    # Check if input file exists
    if not os.path.exists(args.input_fasta):
        print(f"Error: Input file {args.input_fasta} not found", file=sys.stderr)
        sys.exit(1)
    
    # Load the model
    try:
        with open(model_path, 'rb') as fp:
            model = pickle.load(fp)
    except FileNotFoundError:
        print(f"Error: Model file {model_path} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading model: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Load columns reference
    try:
        columns = pd.read_csv(columns_path).iloc[:,1:]
    except FileNotFoundError:
        print(f"Error: Columns file {columns_path} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading columns file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Run BLAST
    output_file = run_blastn(args.input_fasta)
    if output_file is None:
        print("Pathotype-negative")
        sys.exit(0)
    
    # Process BLAST results
    low = False
    try:
        blast_res = pd.read_csv(
            output_file,
            sep="\t",
            names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                   "evalue", "bitscore", "slen"],
        )
        blast_res['cov'] = blast_res['length'] / blast_res['slen']
        blast_res = blast_res[blast_res['cov'] >= 0.6]
        if blast_res.empty:
            low = True
    except pd.errors.EmptyDataError:
        low = True
    except Exception as e:
        print(f"Error processing BLAST results: {e}", file=sys.stderr)
        low = True
    
    if not low:
        # Setup abricate databases first
        if not setup_abricate_databases():
            print("Pathotype-negative")
            cleanup_files(output_file)
            sys.exit(0)
        
        # Run abricate pipeline
        if run_abricate(args.input_fasta):
            try:
                inference_df = pd.read_csv("temp/abr.tsv", sep="\t").iloc[:, 2:]

                # Preprocess the input data
                inference_df = (
                    inference_df.mask(inference_df == ".", 0)
                    .replace({v: 1 for v in inference_df.stack().unique() if v != 0})
                    .mask(inference_df != 0, 1)
                )

                drop_columns = ['EAEC_aggR', 'ETEC_sth', 'STEC_stx1', 'STEC_stx2']
                existing_cols = [c for c in drop_columns if c in inference_df.columns]
                inference_df = inference_df.drop(existing_cols, axis=1)
                inference_df.columns = inference_df.columns.str.replace(r"^ETEC_", "", regex=True)
                inference_df.columns = [col.split('_')[0] for col in inference_df.columns]

                case_map = {}
                new_columns = []
                for original_col in inference_df.columns:
                    lower_col = original_col.lower()
                    if lower_col not in case_map:
                        case_map[lower_col] = original_col
                        new_columns.append(original_col)
                    else:
                        new_columns.append(case_map[lower_col])

                inference_df = (
                    inference_df.set_axis(new_columns, axis='columns')
                    .T.groupby(level=0).sum().T
                    .clip(upper=1)
                    .assign(VF_COUNT=lambda df: df.sum(axis='columns'))
                )

                inference_df = inference_df.reindex(columns=columns.columns, fill_value=0)

                # Drop columns if they exist
                columns_to_drop = ['Sample', 'Pathotype', 'Label']
                existing_columns = [col for col in columns_to_drop if col in inference_df.columns]
                if existing_columns:
                    inference_df = inference_df.drop(existing_columns, axis=1)

                prediction_code = model.predict(inference_df)[0]
                res = "EpiHR" if prediction_code == 0 else "non-EpiHR"

                detected_vfs = [col for col in inference_df.columns if col != 'VF_COUNT' and inference_df[col].values[0] == 1]
                vf_count = int(inference_df['VF_COUNT'].values[0]) if 'VF_COUNT' in inference_df.columns else len(detected_vfs)

                closest_references = []
                if 'blast_res' in locals() and not blast_res.empty:
                    top_matches = blast_res.sort_values(by=['pident', 'cov'], ascending=False).head(3)
                    for _, row in top_matches.iterrows():
                        closest_references.append({
                            "reference_strain": row['sseqid'],
                            "identity_percentage": round(row['pident'], 2),
                            "coverage": round(row['cov'], 2)
                        })

                import json
                output_report = {
                    "input_fasta": args.input_fasta,
                    "prediction": res,
                    "interpretation_basis": {
                        "closest_reference_strains": closest_references,
                        "detected_virulence_genes_count": vf_count,
                        "detected_virulence_genes": detected_vfs
                    }
                }

                print(json.dumps(output_report, indent=4, ensure_ascii=False))

            except Exception as e:
                print(f"Error processing abricate results: {e}", file=sys.stderr)
                import json
                print(json.dumps({"input_fasta": args.input_fasta, "prediction": "Pathotype-negative", "error": str(e)}))
        else:
            import json
            print(json.dumps({"input_fasta": args.input_fasta, "prediction": "Pathotype-negative"}))
    else:
        import json
        print(json.dumps({"input_fasta": args.input_fasta, "prediction": "Pathotype-negative"}))
    
    cleanup_files(output_file)

if __name__ == "__main__":
    main()
