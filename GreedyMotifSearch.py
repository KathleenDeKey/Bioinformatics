# Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
# Input: A string Text, an integer k, and a 4 × k matrix Profile.
# Profile sequence: ACGT
# Output: A Profile-most probable k-mer in Text.

def profile_most_probably_kmer(Dna, k, profile: dict):
    highest_probability = 0
    most_probable_kmer = Dna[:k]
    for i in range(0, len(Dna) - k + 1):
        kmer = Dna[i:i + k]
        kmer_probability = 1
        for position in range(len(kmer)):
            nucleotide = kmer[position]
            kmer_probability *= profile[nucleotide][position]
        if highest_probability < kmer_probability:
            highest_probability = kmer_probability
            most_probable_kmer = kmer
    return most_probable_kmer


# Input: Integers k and t, followed by a space-separated collection of strings Dna.
# Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t). If at any step
# you find more than one Profile-most probable k-mer in a given string, use the one occurring first
from MotifScoring import score
def greedy_motif_search(Dna_list, k, t):
    best_motifs = [Dna[:k] for Dna in Dna_list]
    for i in range(len(Dna_list[0]) - k + 1):
        motif1 = Dna_list[0][i:i + k]
        motifs = [motif1]
        for j in range(1, t):
            profile = create_profile(motifs)
            most_probable_kmer = profile_most_probably_kmer(Dna_list[j], k, profile)
            motifs.append(most_probable_kmer)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


# Helper method to create profile of particular motif in greedy motif search - used Laplace's statictic
from MotifScoring import occurrence
def create_profile(motifs):
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    for col in range(len(motifs[0])):
        nucleotide_count = 0
        for nucleotide in 'ACGT':
            curr_nucleo_occurrence = occurrence(motifs, nucleotide, col)
            curr_nucleo_occurrence += 1
            profile[nucleotide].append(curr_nucleo_occurrence)
            nucleotide_count += curr_nucleo_occurrence
        for nucleotide in 'ACGT':
            profile[nucleotide][col] = profile[nucleotide][col] / nucleotide_count
    return profile

