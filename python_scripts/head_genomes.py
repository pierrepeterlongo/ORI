import pyfastx
import argparse
import sys

if __name__ == "__main__":
    # use argpart to get the file of files
    import argparse
    parser = argparse.ArgumentParser(description='Create a fasta sequence from first N nucleotides of the input fastx file.')
    parser.add_argument('--fastx_file', type=str, help='Path to the file containing paths to sequences.', required=True)
    parser.add_argument('-N', type=int, help='Number of nucleotides to extract from the start of each genome.', default=100)
    args = parser.parse_args()
    
    # open both the file of files and the output file
    with open(args.file_of_files, 'r') as f, open(args.output_file, 'w') as out_f:
        for file in f:
            file = file.strip()
            full_sequence = ""
            try:
                for name, seq in pyfastx.Fasta(file, build_index=False):
                    full_sequence += seq
                if args.N is not None and args.N > 0:
                    full_sequence = full_sequence[:args.N]
                print(f">{name}\n{full_sequence[:args.N]}")
            except Exception as e:
                print(f"Error reading {file}: {e}", file=sys.stderr)