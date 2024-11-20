from Bio import AlignIO
from collections import Counter
import argparse

def get_custom_consensus_from_aln(aln_file, consensus_fa, threshold=0.6):
    try:
        # IUPAC nucleotide codes
        iupac_codes = {
        frozenset(['A']): 'A',
        frozenset(['C']): 'C',
        frozenset(['G']): 'G',
        frozenset(['T']): 'T',
        frozenset(['A', 'G']): 'R',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S',
        frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C']): 'M',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'C', 'G', 'T']): 'N',
    }
        alignment = AlignIO.read(aln_file, "clustal")
        num_sequences = len(alignment) # number of sequences
        alignment_length = alignment.get_alignment_length() # length of each sequences
        consensus = []
        for i in range(alignment_length):
            entire_column = alignment[:, i]
            column = entire_column.replace("-", "") # delete gaps
            counter = Counter(column) # create a speicial dict
            most_common_residue, count = counter.most_common(1)[0]
            if count / num_sequences >= threshold:
                consensus.append(most_common_residue)
            else:
                residue_set = frozenset([residue for residue in counter if residue != '-']) # delele '-' from column
                iupac_code = iupac_codes.get(residue_set, 'N')
                consensus.append(iupac_code)
        consensus_seq = ''.join(consensus) # make a string
        consensus_name = aln_file.split("/")[-1][:-4]
        with open(consensus_fa, "w") as f: # save in fasta format
            f.write(f'>{consensus_name}_consensus' + "\n")
            f.write(consensus_seq)
    except ValueError as e:
        print(f"Error occured with file processing {aln_file}: {e}")
    except FileNotFoundError as e:
        print(f"Error occured with file processing {consensus_fa}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a consensus sequence from an alignment file.")
    # Define command line arguments
    parser.add_argument("-aln", type=str, help="Input alignment file in Clustal format")
    parser.add_argument("-cons", help="Output file for the consensus sequence in FASTA format")
    parser.add_argument("-t", type=float, default=0.6, help="Threshold for consensus (default: 0.6)")
    
    # Parse the arguments
    args = parser.parse_args()

    get_custom_consensus_from_aln(args.aln, args.cons, args.t)
