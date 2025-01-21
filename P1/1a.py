class Read_map():
    def __init__(self, reference, reads, threshold, paired = False):
        self.reference = reference
        self.reads = reads
        self.threshold = threshold
        self.paired = paired

        if paired:
            self.paired_dist = self.get_paired_dist()


    def get_paired_dist(self):
        # calculate the distance between the paired reads
        # assume the distance is same for all paired reads

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
        dup_list = {}
        for read in self.reads:
            start, possible_mutation_loc = self.align_read(read)

            # if the read is aligned
            if start > 0:
                for j in possible_mutation_loc:
                    if start + j not in dup_list:
                        dup_list[start + j] = [read[j]]
                    else:
                        dup_list[start + j].append(read[j])

        for pos, value in dup_list.items():
            # take the most common mutation in value
            mutation = max(value, key = value.count)
            mutation_list.write(
                '>S' + str(pos) + ' ' + self.reference[pos] + ' ' + mutation + '\n')
                        
        mutation_list.close()

    def create_position_map(self):
        pos_map = {}
        len_frag = self.len_read//5

        for read in self.reads:
            frag_list = [read[0:len_frag], 
                        read[len_frag:2*len_frag], 
                        read[2*len_frag:3*len_frag], 
                        read[3*len_frag:4*len_frag], 
                        read[4*len_frag:5*len_frag]]
            
            # for each fragment, find the position in the reference and store it in pos_map
            for i, frag in enumerate(frag_list):
                # search the reference for the fragment
                search_ref = self.reference
                while search_ref != '':
                    pos = search_ref.find(frag)
                    if pos != -1:
                        search_ref = search_ref[pos+len_frag:]
                        if pos in pos_map:
                            pos_map[frag].append((pos, i))
                        else:
                            pos_map[frag] = [(pos, i)]
                    else:
                        break
        return pos_map
    
    def read_map_indexing(self):
        # create predcited_mutations.txt
        # use the position map to find the mutations
        mutation_list = open('predicted_mutations.txt', 'w')
        pos_map = self.create_position_map()
        dup_list = set()
        for frag, (loc, idx) in pos_map.items():
            pass
            

            

            




def main():
    reference = open('project1a_reference_genome.fasta', 'r')
    reference = ''.join(reference.readlines()[1:]).replace('\n', '')

    reads = open('project1a_with_error_paired_reads.fasta', 'r')
    reads = [line.replace('\n', '') for line in reads.readlines() if line[0] != '>']

    # read_map
    mapping = Read_map(reference, reads, 5, True)
    mapping.read_map_sliding_window()

    # mapping.read_map_sliding_window()


if __name__ == '__main__':
    main()
