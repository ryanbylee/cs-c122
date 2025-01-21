def pattern_count(text: str, pattern: str) -> int:
    '''
    Input: Strings Text and Pattern.
    Output: Count(Text, Pattern).
    '''

    count = 0
    
    for i in range(len(text) - len(pattern) + 1):
        if text[i : i + len(pattern)] == pattern:
            count += 1

    return count

def MaxMap(map):
    max = 0
    max_res = []
    for val in map.values():
        if val >= max:
            max = val
    
    for key, val in map.items():
        if val == max:
            max_res.append(key)
    return max_res

def frequent_words(text: str, k: int) -> list[str]:
    """Find the most frequent k-mers in a given text."""
    freq = {}

    i = 0
    while i + k - 1 < len(text):
        kmer = text[i: i + k]
        freq[kmer] = freq.get(kmer, 0) + 1

        i += 1

    
    return MaxMap(freq)


def reverse_complement(pattern: str) -> str:
    """Calculate the reverse complement of a DNA pattern."""
    comp_map = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

    new_str = ''
    for char in pattern:
        new_str += (comp_map[char])

    return new_str[::-1] 


def pattern_matching(pattern: str, genome: str) -> list[int]:
    """Find all occurrences of a pattern in a genome."""
    len_pattern = len(pattern)
    res = []
    i = 0
    while i < len(genome) - len_pattern + 1:
        if genome[i: i + len_pattern] == pattern:
            res.append(i)
        i += 1
    return res
    
def find_clumps(genome: str, k: int, l: int, t: int) -> list[str]:
    """Find patterns forming clumps in a genome."""
    res = []
    i = 0
    while i < len(genome) - l + 1:
        window = genome[i: i + l]
        freq = {}
        j = 0
        while j + k - 1 < len(window):
            kmer = window[j: j + k]
            freq[kmer] = freq.get(kmer, 0) + 1

            j += 1

        for key, val in freq.items():
            if val >= t:
                res.append(key)
        i += 1
    return set(res)


def minimum_skew(genome: str) -> list[int]:
    """Find positions in a genome where the skew diagram attains a minimum."""
    # g inc, c dec
    res = []
    nuc_to_num = {
        'G': 1,
        'C': -1
    }
    skew_track = []
    cur_val = 0
    min = 0
    for i in range(len(genome)):
        if genome[i] in nuc_to_num:
            cur_val += nuc_to_num[genome[i]]
            
        skew_track.append(cur_val)
        if cur_val < min:
            min = cur_val
    # traverse thru skew_track and locate min
    for i in range(1, len(skew_track)):
        if skew_track[i] == min:
            res.append(i + 1)
    return res


        
def hamming_distance(p: str, q: str) -> int:
    """Calculate the Hamming distance between two strings."""
    count = 0
    for i, j in zip(p, q):
        if i != j:
            count += 1
    return count
    
def approximate_pattern_matching(pattern: str, text: str, d: int) -> list[int]:
    """Find all starting positions where Pattern appears as a substring of Text with at most d mismatches."""
    len_pattern = len(pattern)
    res = []
    i = 0
    while i < len(text) - len_pattern + 1:
        if hamming_distance(text[i: i + len_pattern],  pattern) <= d:
            res.append(i)
        i += 1
    return res


def approximate_pattern_count(text: str, pattern: str, d: int) -> int:
    """Count the occurrences of a pattern in a text, allowing for up to d mismatches."""
    count = 0
    
    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(text[i : i + len(pattern)], pattern) <= d:
            count += 1

    return count

def frequent_words_with_mismatches(text: str, k: int, d: int) -> list[str]:
    """Find the most frequent k-mers with up to d mismatches in a text."""
    freq = {}
    n = len(text)

    i = 0
    while i + k - 1 < n:
        pattern = text[i: i+k]
        neighborhood = neighbors(pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            freq[neighbor] = freq.get(neighbor, 0) + 1
        i += 1

    return MaxMap(freq)




def neighbors(s: str, d: int) -> list[str]:
    """Generate neighbors of a string within a given Hamming distance."""

    nucs = ['A', 'C', 'G', 'T']
    suffix = s[1:]

    if d == 0:
        return [s]
    if len(s) == 1:
        return nucs
    
    neighborhood = set()
    suffixNeighbors = neighbors(suffix, d)

    for sufNeighbor in suffixNeighbors:
        if hamming_distance(suffix, sufNeighbor) < d:
            for nuc in nucs:
                neighborhood.add(nuc + sufNeighbor)
        else:
            neighborhood.add(s[0] + sufNeighbor)
    return list(neighborhood)

def main():
    print(frequent_words_with_mismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1))
    

if __name__ == '__main__':
    main()