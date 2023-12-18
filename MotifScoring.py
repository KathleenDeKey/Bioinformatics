import math


# Use entropy equation to score a list of Dna strings' variations
def score(motifs):
    motifs = [[nucleo for nucleo in motif] for motif in motifs]
    counts = {'A': [], 'T': [], 'C': [], 'G': []}
    for col in range(len(motifs[0])):
        for nucleotide in 'ATCG':
            counts[nucleotide].append(occurrence(motifs, nucleotide, col))
    total_entropy = 0
    for i in range(len(motifs[0])):
        entropy = 0
        for nucleotide in 'ATCG':
            occurence = counts[nucleotide][i]
            if occurence != 0:
                entropy += occurence * math.log2(occurence)
        entropy = -entropy
        total_entropy += entropy
    return total_entropy


# Calculate the occurrence of a nucleotide in the specific position
def occurrence(motifs, nucleotide, col):
    occurrence = 0
    for row in motifs:
        curr_nucleo = row[col]
        if curr_nucleo.upper() == nucleotide:
            occurrence += 1
    return occurrence

