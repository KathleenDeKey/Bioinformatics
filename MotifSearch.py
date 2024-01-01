# Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
# Input: A collection of strings Dna, and integers k and d.
# Output: All (k, d)-motifs in Dna.
def motif_enumeration(Dna_list, k, d):
    patterns = set()
    for Dna_sequence in Dna_list:
        for i in range(0, len(Dna_sequence) - k + 1):
            k_mer = Dna_sequence[i:i + k]
            neighbors = generate_neighbors(k_mer, d)
            for neighbor in neighbors:
                if is_motif_in_Dna(neighbor, Dna_list, d):
                    patterns.add(neighbor)
        return patterns


# Finds all k_mers that differ from 'pattern' by at most 'd' mismatches
def generate_neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A','C','G','T'}
    neighbors = set()
    suffix_neighbors = generate_neighbors(pattern[1:], d)
    for neighbor in suffix_neighbors:
        if hamming_distance(pattern[1:], neighbor) < d:
            for nucleotide in 'ACGT':
                neighbors.add(nucleotide + neighbor)
        else:
            neighbors.add(pattern[0] + neighbor)
    return neighbors


# Calculate the hamming distance between two sequences of equal length
def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


# Find if a pattern appears in all Dna strings with at most d mismatches
def is_motif_in_Dna(k_mer, Dna_list, d):
    k = len(k_mer)
    for Dna_sequence in Dna_list:
        found = False
        for i in range(len(Dna_sequence) - k + 1):
            temp_k_mer = Dna_sequence[i:i+k]
            if hamming_distance(temp_k_mer, k_mer) <= d:
                found = True
                break
        if not found:
            return False
    return True


