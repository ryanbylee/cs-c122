# TODO:
# create kmers of each read and them in some data structure?
# create hash function for bloom filter
# create bloom filter for each genome


from tqdm import tqdm
from bitarray import bitarray

'''
useful info:
l = num of hash function
m = length of filter
n = num of elements to be added

optimal l = 6.64 ~ 7
optimal m/n (ratio of length vs. num of reads) = 9.56 ~ 10

30 kmers per read

4mil length for the filter
'''


class Read_mapper():
    def __init__(self, genomes, reads):
        self.genomes = self.process_genomes(genomes)
        self.hash_functions = self.generate_hash_functions(7)
        self.bloom_filters = []

        print('\ncreating bloom filters...')
        for genome in tqdm(self.genomes.values()):
            self.bloom_filters.append(self.create_bloom_filter(genome, self.hash_functions))

        self.reads = reads

    def process_genomes(self, genomes):
        # process genomes
        # create dictionary of genome number to genome sequence

        genome_dict = {}
        print('\nprocessing genomes...')
        for genome in tqdm(genomes):
            genome_dict[genome[0]] = genome[1]

        return genome_dict
    

    def kmerize(self, sequence):
        # create kmers of the read
        # return a list of kmers

        kmer_list = []
        for i in range(len(sequence) - 15 + 1):
            window = sequence[i:i+15]
            if window not in kmer_list:
                kmer_list.append(sequence[i:i+15])

        return kmer_list


    def generate_hash_functions(self, num_hashes):
        # create {num_hashes} amount of universal hash functions
        # return a list of hash functions
        print('\ngenerating hash functions...')
        # h_a,b = ((ax + b) % p) % m
        hash_functions = []
        a, b = 1, 1
        for _ in tqdm(range(num_hashes)):
            hash_functions.append(lambda x: ((a * hash(x) + b) % 400009) % 4000000)
            a += 1
            b += 1  

        return hash_functions
    
    def create_bloom_filter(self, genome, hash_functions):

        bloom_filter = bitarray(4000000)

        kmers = self.kmerize(genome)

        for kmer in kmers:
            for hash_function in hash_functions:
                bloom_filter[hash_function(kmer)] = 1

        return bloom_filter
    
    def check_bloom_filter(self, read, bloom_filter):

        kmers = self.kmerize(read)
        count = 0
        for kmer in kmers:
            match_all_hash = True
            for hash_function in self.hash_functions:
                if bloom_filter[hash_function(kmer)] == 0:
                    match_all_hash = False
            if match_all_hash:
                count += 1
        return count
    
    def map_reads(self):

        read_map = {}
        for i, read in enumerate(tqdm(self.reads)):
            for j, bloom_filter in enumerate(self.bloom_filters):

                count = self.check_bloom_filter(read, bloom_filter)
                if count > 0:
                    if i in read_map:
                        read_map[i].append((j, count))  
                    else:
                        read_map[i] = [(j, count)]

        return self.take_consensus(read_map)
    

    def create_freq_map(self, read_map):
        # helper function for 1d
        freq = {}
        for read, genome_src in read_map.items():
            freq[genome_src] = freq.get(genome_src, 0) + 1
        
        return freq
    
    def take_consensus(self, read_map):
        
        print('\ntaking consensus...')
        consensus_map = {}
        for read, genome_list in tqdm(read_map.items()):
            max_genome = -1
            max_count = -1
            for genome, count in genome_list:
                if count > max_count:
                    max_count = count
                    max_genome = genome
            consensus_map[read] = max_genome

        return consensus_map
    
    def output_res(self, consensus_map):

        res = open('project1c-files/answers.txt', 'w')
        print('\nwriting results to file...')
        for read, genome in tqdm(consensus_map.items()):
            res.write(f'>read_{read}\tGenome_Number_{genome}\n')
        res.close()
    
def main():

    genomes = []
    for i in range(0, 1000):
        genome = open(f'project1d-files/project1d_genome_{i}.fasta', 'r')
        genome = ''.join(genome.readlines()[1:]).replace('\n', '')

        genomes.append((i, genome))


    reads = open('project1d-files/project1d_reads.fasta', 'r')
    reads = [line.replace('\n', '') for line in reads.readlines() if line[0] != '>']
    # for 1d, take only 10000 reads
    reads = reads[0:100000]



    mapper = Read_mapper(genomes, reads)
    read_map = mapper.map_reads()

    freq = mapper.create_freq_map(read_map)

    print(freq.keys())
    # mapper.output_res(read_map)



if __name__ == '__main__':    
    main()