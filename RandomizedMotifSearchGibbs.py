# Gibbs sampler that randomly selects a string of DNA to ignore
# Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
# Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts.
# Remember to use pseudocounts!
from RandomizedMotifSearch import generate_random_kmer, score
import random
from GreedyMotifSearch import create_profile
def gibbs_sampler(Dna_list, k, t, N):
    random_motifs = generate_random_kmer(Dna_list, k)
    best_motifs = random_motifs
    j = 1
    while j < N:
        i = random.randint(0, t - 1)
        temp_motifs = [motif for motif in best_motifs]
        profile = create_profile(temp_motifs[:i] + temp_motifs[i + 1:])
        motifi = profile_randomly_generated_kmer(Dna_list[i], k, profile)
        temp_motifs[i] = motifi
        if score(temp_motifs, k, t) < score(best_motifs, k, t):
            # print('temp motif: ', temp_motifs)
            best_motifs = temp_motifs
            # print('best motif: ', best_motifs)
        j += 1
        # print(j)
    return best_motifs


# Random number generator Random(p1, …, pn) that returns an index between 0 and n-1 according to the probability
def random_generator(probabilities):
    total_prob = sum(probabilities)
    if total_prob <= 0:
        raise ValueError("Probabilities should sum to a positive constant.")
    normalized_probs = [p / total_prob for p in probabilities]
    rand_num = random.uniform(0, 1)
    cumulative_prob = 0
    for i, prob in enumerate(normalized_probs):
        cumulative_prob += prob
        if rand_num <= cumulative_prob:
            return i


# Returns a kmer that is randomly selected by using the probability distribution
def profile_randomly_generated_kmer(Dna, k, profile):
    probabilities = []
    for i in range(0, len(Dna) - k + 1):
        kmer = Dna[i:i + k]
        kmer_probability = 1
        for position in range(len(kmer)):
            nucleotide = kmer[position]
            kmer_probability *= profile[nucleotide][position]
        probabilities.append(kmer_probability)
    # print("probabilities: ", probabilities)
    rand_kmer_index = random_generator(probabilities)
    # print("rand kmer index:", rand_kmer_index)
    rand_kmer = Dna[rand_kmer_index:rand_kmer_index + k]
    # print(rand_kmer)
    return rand_kmer


# Run the randomized motif search with gibbs a given number of times to get an approximation
def result(Dna, k, t, N):
    best_score = float("inf")
    best_motifs = []
    i = 0
    while i < 35:
        curr_best_motifs = list(gibbs_sampler(Dna, k, t, N))
        curr_score = score(curr_best_motifs, k, t)
        if curr_score < best_score:
            # print('in if loop')
            best_score = curr_score
            best_motifs = curr_best_motifs
            i = 0
        i += 1
    return best_motifs

