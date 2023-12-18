# The first potential issue with implementing MedianString is writing a function to compute
# d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai), the sum of distances between Pattern and each string in
# Dna = {Dna1, ..., Dnat}.
# Input: A string Pattern followed by a collection of space-separated strings Dna.
# Output: d(Pattern, Dna)
from MotifSearch import hamming_distance


def distance_between_pattern_and_strings(pattern, dna_list):
    k = len(pattern)
    distance = 0
    for dna in dna_list:
        min_hamming_distance = float("inf")
        for i in range(0, len(dna) - k + 1):
            k_mer = dna[i:i + k]
            temp_hamming = hamming_distance(pattern, k_mer)
            min_hamming_distance = min(min_hamming_distance, temp_hamming)
        distance += min_hamming_distance
    return distance


# Input: An integer k, followed by a space-separated collection of strings Dna.
# Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
# (If there are multiple such strings Pattern, then you may return any one.)
def median_string(dna_list, k):
    distance = float('inf')
    patterns = all_strings(k)
    median = ""
    for pattern in patterns:
        curr_distance = distance_between_pattern_and_strings(pattern, dna_list)
        if distance > curr_distance:
            distance = curr_distance
            median = pattern
    return median


# Returns an array containing all strings of length k.
from MotifSearch import generate_neighbors


def all_strings(k):
    pattern = 'A' * k
    return generate_neighbors(pattern, k)


