class Read_map():
    def __init__(self, reference, read, threshold):
        self.reference = reference
        self.reads = read
        self.threshold = threshold
        self.len_read = 50

    def read_map_sliding_window(self):
        # sliding window
        # create predcited_mutations.txt

        mutation_list = open('predicted_mutations.txt', 'w')
        dup_list = {}
        for read in self.reads:
            
            read_length = len(read)
            if read_length < self.len_read:
                continue
            for i in range(len(self.reference) - read_length + 1):
                window = self.reference[i:i+read_length]
                count = 0
                for j in range(read_length):
                    if read[j] != window[j]:
                        count += 1

                if count <= self.threshold:
                    for j in range(read_length):

                        if read[j] != window[j]:
                            if i + j not in dup_list:
                                dup_list[i + j] = [read[j]]
                            else:
                                dup_list[i + j].append(read[j])
                        
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
    reference = open('sample_reference_genome.fasta', 'r')
    reference = ''.join(reference.readlines()[1:]).replace('\n', '')

    sample_no_error_single = open('sample_no_error_paired_reads.fasta', 'r')
    sample_no_error_single = sample_no_error_single.readlines()
    sample_no_error_single = [line.replace('\n', '') for line in sample_no_error_single if line[0] != '>']

    # read_map
    mapping = Read_map(reference, sample_no_error_single, 3)
    mapping.read_map_sliding_window()

    # mapping.read_map_sliding_window()


    # sanity check
    pred = open('predicted_mutations.txt', 'r')
    ans = open('sample_mutations.txt', 'r')

    pred = pred.readlines()
    ans = ans.readlines()

    ans_filtered = [line for line in ans if line[0:2] == '>S']

    for line in ans_filtered:
        if line not in pred:
            print("not matched", line)


if __name__ == '__main__':
    main()
