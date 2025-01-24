import itertools

class Read_mapper():
    def __init__(self, reference, reads, threshold, paired = False):
        self.reference = reference
        self.reads = reads
        self.threshold = threshold
        self.paired = paired

        # if paired:
        #     self.paired_dist = self.get_paired_dist()


    def get_paired_dist(self):
        # calculate the distance between the paired reads
        # assume the distance is same for all paired reads
        # useful info, but not necessary

        first_read = self.reads[6]
        second_read = self.reads[7]

        pos1 = self.align_read(first_read)
        pos2 = self.align_read(second_read)

        distance = pos2 - (pos1 + len(first_read))
        return distance
    

    def align_read(self, read):
        # align the reads to the reference genome
        # return the position of the read in the reference genome
        read_length = len(read)
        
        # ignore reads that are not of length 50
        if read_length != 50:
            return -1, None
        
        for i in range(len(self.reference) - read_length + 1):
            window = self.reference[i:i+read_length]
            count = 0
            possible_mutation_loc = []
            for j in range(read_length):
                if read[j] != window[j]:
                    count += 1
                    possible_mutation_loc.append(j)

            if count <= self.threshold:
                return i, possible_mutation_loc
            
        return -1, None
    
    def read_map_sliding_window(self):
        # sliding window
        # create predcited_mutations.txt
        mutation_list = open('predicted_mutations.txt', 'w')
        overlapping_mutations = {}
        for read in self.reads:
            start, possible_mutation_loc = self.align_read(read)

            # if the read is aligned
            if start > 0:
                for j in possible_mutation_loc:
                    if start + j not in overlapping_mutations:
                        overlapping_mutations[start + j] = [read[j]]
                    else:
                        overlapping_mutations[start + j].append(read[j])

        for pos, value in overlapping_mutations.items():
            # if the mutation is not present in at least 3 reads, ignore it
            if len(value) < 3:
                continue

            # take the majority mutation in value
            mutation = max(value, key = value.count)

            mutation_list.write(
                '>S' + str(pos) + ' ' + self.reference[pos] + ' ' + mutation + '\n')
                        
        mutation_list.close()

    # not used; meant to be a helper function for hashmap approach
    def create_position_map(self):
        pos_map = {}
        len_frag = 10

        # create the keys for the dictionary, which are all possible strings of length 10 that consist A, T, C, G (~4^10 entries)
        for string in itertools.product('AGCT', repeat=10):
            pos_map[''.join(string)] = []

        
        for string in pos_map.keys():
            # locate the position of the fragment in the reference genome and populate the values of the dictionary
            for i in range(len(self.reference) - len_frag + 1):
                if self.reference[i:i+len_frag] == string:
                    pos_map[string].append(i)

        for read in self.reads:
            # split the read into 5 fragments
            frag_list = [read[0:len_frag], 
                        read[len_frag:2*len_frag], 
                        read[2*len_frag:3*len_frag], 
                        read[3*len_frag:4*len_frag], 
                        read[4*len_frag:5*len_frag]]
            
            for i, frag in enumerate(frag_list):
                # search for the fragment in pos_map
                if frag in pos_map:
                    print(frag + ' found at: ' + str(pos_map[frag]))

        return pos_map



def main():
    reference = open('project1a_reference_genome.fasta', 'r')
    reference = ''.join(reference.readlines()[1:]).replace('\n', '')

    reads = open('project1a_with_error_paired_reads.fasta', 'r')
    reads = [line.replace('\n', '') for line in reads.readlines() if line[0] != '>']

    # read_map
    mapping = Read_mapper(reference, reads, 4, True)
    mapping.read_map_sliding_window()
    # mapping.create_position_map()



if __name__ == '__main__':
    main()
