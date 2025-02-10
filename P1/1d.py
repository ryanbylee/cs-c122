import itertools
from tqdm import tqdm


'''
1d approach: use minimizer and needleman-wunsch

'''

class Metagenome():
    def __init__(self, genomes, reads, kmer_size = 21, minimizer_size = 15):
        self.genomes = genomes
        self.kmer_size = kmer_size
        self.minimizer_size = minimizer_size
        # minimizer_maps: list of minimizer maps for each genome, key: minimizer, value: list of positions
        self.minimizer_maps = [self.create_minimizer_map(genome) for genome in tqdm(genomes, desc='creating minimizer maps...')]
        self.reads = reads

    

    def calculate_minimizer(self, seq):
        # calculate the minimizer of a sequence

        minimizer = min([seq[i:i+self.minimizer_size] for i in range(len(seq) - self.minimizer_size + 1)])

        return minimizer
    
    def needleman_wunsch(self, read, ref):
        # Needleman-Wunsch Algorithm, global alignment using dynamic programming

        read_length = len(read)

        # initialize dp table
        dp = [ [0 for i in range(read_length + 1)] for j in range(len(ref) + 1)]

        for i in range(1, read_length + 1):
            dp[i][0] = i
        for j in range(1, len(ref) + 1):
            dp[0][j] = j
        # list of tuples in form of (type of operation, position)

        def backtrack(dp, read, i, j):
            # backtrack the dp table to find the alignment

            # base case
            if i == 0 and j == 0:
                return
                                
            insert = dp[i][j-1] + 1
            delete = dp[i-1][j] + 1
            match = dp[i-1][j-1] + (0 if read[i-1] == ref[j-1] else 1)

            best_score = min(insert, delete, match)

    
            if best_score == insert:
                return backtrack(dp, read, i, j-1)
            elif best_score == delete:
                return backtrack(dp, read, i-1, j)
            elif best_score == match:
                return backtrack(dp, read, i-1, j-1)

            

        # fill the dp table
        for i in range(1, read_length + 1):
            for j in range(1, len(ref) + 1):
                match = dp[i-1][j-1] + (0 if read[i-1] == ref[j-1] else 1)
                delete = dp[i-1][j] + 1
                insert = dp[i][j-1] + 1
                dp[i][j] = min(match, delete, insert)

        backtrack(dp, read, read_length, len(ref))

        return dp[read_length][len(ref)]
            

    def create_minimizer_map(self, genome):

        minimizer_map = {}
       
        for i in range(len(genome) - self.kmer_size + 1):
            kmer = genome[i:i+self.kmer_size]

            # alphabetical minimum
            minimizer = self.calculate_minimizer(kmer)

            if minimizer not in minimizer_map:
                minimizer_map[minimizer] = [i]
            else:
                minimizer_map[minimizer].append(i)

        return minimizer_map

    
    def find_candidate_pos(self, read, minimizer_map):

        # make kmer of the read, and calculate minimizer for each kmer
        read_length = len(read)
        kmers = [read[i:i+self.kmer_size] for i in range(read_length - self.kmer_size + 1)]
        minimizers = list(set([self.calculate_minimizer(kmer) for kmer in kmers]))

        candidate_pos = []

        for minimizer in minimizers:
            if minimizer in minimizer_map:
                candidate_pos.extend(minimizer_map[minimizer])


        return candidate_pos
    
    def count_occurance(self, read):
        # for read, go thru all minimizer map, get candidate positions, and log the count of positions that require < 7 mutations
        # return the genome number with the most occurances
        max_count = 0
        max_genome = -1
        for i, minimizer_map in enumerate(self.minimizer_maps):
            candidate_pos = self.find_candidate_pos(read, minimizer_map)
            count = 0

            for pos in candidate_pos:
                ref_window = self.genomes[i][pos: pos + len(read)]
                if len(ref_window) != len(read):
                    continue
                score = self.needleman_wunsch(read, ref_window)

                if score < 7:
                    count += 1

            if count > max_count:
                max_count = count
                max_genome = i

        return max_genome

def main():
    genomes = []

    # take first 100 genomes
    for i in range(0, 10):
        genome = open(f'project1c_sample/project1c_sample_genome_{i}.fasta', 'r')
        genome = ''.join(genome.readlines()[1:]).replace('\n', '')

        genomes.append(genome)

    reads = open('project1c_sample/project1c_sample_reads.fasta', 'r')
    reads = [line.replace('\n', '') for line in tqdm(reads.readlines(), desc='loading reads...') if line[0] != '>']

    # take only the first 2000 reads

    # read_map
    mapping = Metagenome(genomes, reads)


    res = open('project1c_sample/answers.txt', 'w')
    i = 0
    for read in tqdm(reads, desc='mapping reads...'):
        res.write(f'>read_{i}\tGenome_Number_{mapping.count_occurance(read)}\n')
        i += 1
    res.close()
    



if __name__ == '__main__':
    main()
