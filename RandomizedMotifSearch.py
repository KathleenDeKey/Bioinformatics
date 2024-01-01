# Input: A collection of strings Dna, and integers k and t.
# Output: A collection of strings resulting from running RANDOMIZEDMOTIFSEARCH(Dna, k, t) 1000
# times. Remember to use pseudocounts
from GreedyMotifSearch import create_profile
def randomized_motif_search(Dna_list, k, t):
    random_motifs = generate_random_kmer(Dna_list, k)
    best_motifs = random_motifs
    while True:
        profile = create_profile(best_motifs)
        probable_motifs = generate_probable_motifs(profile, Dna_list, k)
        if score(probable_motifs, k, t) < score(best_motifs, k, t):
            best_motifs = probable_motifs
        else:
            return best_motifs

# Helper methods to generate random kmers
import random
def random_kmer(dna, k):
    start_index = random.randint(0, len(dna) - k)
    return dna[start_index:start_index + k]


def generate_random_kmer(Dna_list, k):
    random_kmers = [random_kmer(dna, k) for dna in Dna_list]
    return random_kmers


# Helper method to find the list of most probable kmers as motifs
from GreedyMotifSearch import profile_most_probably_kmer
def generate_probable_motifs(profile, Dna_list, k):
    motifs = []
    for dna in Dna_list:
        probable_kmer = profile_most_probably_kmer(dna, k, profile)
        motifs.append(probable_kmer)
    return motifs


# Helper Method to Score Motifs based on the deviation simply using counts
import numpy as np
from collections import Counter

def score(motifs, k, t):
    List = np.asarray([list(motif) for motif in motifs])
    score = 0
    for i in range(k):
        c = Counter(List[:, i])
        score += t - c.most_common(1)[0][1]
    return score


# Run the randomized motif search a given number of times to get an approximation
def result(k, t, Dna):
    best_score = float("inf")
    best_motifs = []
    i = 0
    while i < 2000:
        curr_best_motifs = list(randomized_motif_search(Dna, k, t))
        curr_score = score(curr_best_motifs, k, t)
        if curr_score < best_score:
            best_score = curr_score
            best_motifs = curr_best_motifs
            i = 0
        i += 1
    return best_motifs


