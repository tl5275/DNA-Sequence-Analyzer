import numpy as np

# Genetic Code Dictionary for optional translation
genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def count_base_pairs(sequence):
    sequence = sequence.upper()
    return {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }

def identify_mutations(seq1, seq2):
    seq1, seq2 = seq1.upper(), seq2.upper()
    min_len = min(len(seq1), len(seq2))
    arr1 = np.array(list(seq1[:min_len]))
    arr2 = np.array(list(seq2[:min_len]))
    diffs = arr1 != arr2
    return [(i, arr1[i], arr2[i]) for i in np.where(diffs)[0]]

def translate_dna(dna_seq):
    dna_seq = dna_seq.upper()
    protein = ''
    for i in range(0, len(dna_seq)-2, 3):
        codon = dna_seq[i:i+3]
        protein += genetic_code.get(codon, 'X')  # 'X' for unknown/invalid codon
    return protein

# Main interactive block
if __name__ == "__main__":
    print("=== NA Sequence Analyzer ===")
    seq1 = input("Enter Sequence 1 (DNA): ").strip().upper()
    seq2 = input("Enter Sequence 2 (DNA): ").strip().upper()

    print("\n--- Base Pair Counts for Sequence 1 ---")
    base_counts = count_base_pairs(seq1)
    for base, count in base_counts.items():
        print(f"{base}: {count}")

    print("\n--- Mutation Locations ---")
    mutations = identify_mutations(seq1, seq2)
    if mutations:
        for idx, base1, base2 in mutations:
            print(f"Position {idx}: {base1} -> {base2}")
    else:
        print("No mutations found.")

    print("\n--- DNA to Amino Acid Translation (Sequence 1) ---")
    protein = translate_dna(seq1)
    print("Protein:", protein)
