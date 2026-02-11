import pyfastx
import argparse
import sys

if __name__ == "__main__":
    # use argpart to get the file of files
    import argparse
    parser = argparse.ArgumentParser(description='Create a fasta sequence from first N nucleotides of the input fastx file.')
    parser.add_argument('--fastx_file', type=str, help='Path to the file containing paths to sequences.', required=True)
    parser.add_argument('-N', type=int, help='Number of nucleotides to extract from the start of each genome.', default=100)
    parser.add_argument('--start_from', type=int, help='Where to start in each genome.', default=0)
    args = parser.parse_args()
    
    full_sequence = ""
    try:
        for name, seq in pyfastx.Fasta(args.fastx_file, build_index=False):
            full_sequence += seq
            if len(full_sequence) > args.start_from + args.N: break
        print(f">{name}\n{full_sequence[args.start_from:args.start_from+args.N]}")
    except Exception as e:
        print(f"Error reading {args.fastx_file}: {e}", file=sys.stderr)