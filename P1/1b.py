import itertools
from tqdm import tqdm

class Read_mapper():
    def __init__(self, reference, reads, threshold, approach):
        self.reference = reference
        self.reads = reads
        self.threshold = threshold
        self.approach = approach


    def get_paired_dist(self):
        # calculate the distance between the paired reads
        # assume the distance is same for all paired reads
        # useful info, but not necessary

        first_read = self.reads[6]
        second_read = self.reads[7]

        pos1 = self.slide_window(first_read)
        pos2 = self.slide_window(second_read)

        distance = pos2 - (pos1 + len(first_read))
        return distance
    

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
        possible_indel_sub_loc = []

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
                possible_indel_sub_loc.append(('>D', j-1, '_', ref[j-1]))
                return backtrack(dp, read, i, j-1)


            elif best_score == delete:
                possible_indel_sub_loc.append(('>I', i-1, '_', read[i-1]))
                return backtrack(dp, read, i-1, j)
            elif best_score == match:
                # if the characters are the same, do nothing
                if read[i-1] != ref[j-1]:
                    possible_indel_sub_loc.append(('>S', i-1, ref[j-1], read[i-1]))
                return backtrack(dp, read, i-1, j-1)

            

        # fill the dp table
        for i in range(1, read_length + 1):
            for j in range(1, len(ref) + 1):
                match = dp[i-1][j-1] + (0 if read[i-1] == ref[j-1] else 1)
                delete = dp[i-1][j] + 1
                insert = dp[i][j-1] + 1
                dp[i][j] = min(match, delete, insert)

        backtrack(dp, read, read_length, len(ref))

        return dp[read_length][len(ref)], possible_indel_sub_loc[::-1]
            




    def slide_window(self, read):
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
                    possible_mutation_loc.append(('>S', j, window[j], read[j]))

            if count <= self.threshold:
                return i, possible_mutation_loc
            
        return -1, None
    
    def read_map(self):
        # sliding window
        # create predcited_mutations.txt
        mutation_list = open('predicted_mutations.txt', 'w')
        overlapping_mutations = {}
        pos_map = self.create_position_map()
        for read in tqdm(self.reads):

            if self.approach == 'sliding_window':
                # sliding window approach
                start, possible_mutation_loc = self.slide_window(read)

            elif self.approach == 'dp':
                # dynamic programming approach
                
                start, possible_mutation_loc = self.find_best_pos(read, pos_map)
            else:
                raise NotImplementedError(f'approach {self.approach} not implemented')
            
            # if the read is aligned
            if start > 0:
                for op, loc, ref_char, char in possible_mutation_loc:
                    if start + loc not in overlapping_mutations:
                        overlapping_mutations[start + loc] = [(op, start + loc, ref_char, char)]
                    else:
                        overlapping_mutations[start + loc].append((op, start + loc, ref_char, char))

        for pos, value in overlapping_mutations.items():
            # if the mutation is not present in at least 3 reads, ignore it
            if len(value) < 3:
                continue

            # take the majority mutation in value
            mutation = max(value, key = lambda x: x[3])
            
            ref_char = ' ' + mutation[2] + ' ' if mutation[0] == '>S' else ' '

            mutation_list.write(
                mutation[0] + str(pos) + ref_char + mutation[3] + '\n')
                        
        mutation_list.close()

    def create_position_map(self):
        # create a dictionary of all 15-length fragments in the reference genome and their positions
        pos_map = {}
        len_frag = 15

        for i in range(len(self.reference) - len_frag + 1):
            frag = self.reference[i:i+len_frag]
            if frag not in pos_map:
                pos_map[frag] = [i]
            else:
                pos_map[frag].append(i)

        return pos_map

    def find_possible_pos(self, ref_pos_map, read):
        # locate the read in the reference genome

        # divide the read into 15-length fragments and discard the last fragment if it is less than 15
        len_frag = 15
        read_frags = [read[i:i+len_frag] for i in range(0, len(read), len_frag)]

        candidate_pos = [] # store starting positions and the corresponding fragments from reference genome
        # since each read is divided into 3 fragments and the error threshold is 2, there must be at least one fragment out of 3 that matches the reference genome
        for i, frag in enumerate(read_frags):
            if frag in ref_pos_map:
                # 1. the first fragment matches the reference genome
                if i == 0:
                    for pos in ref_pos_map[frag]:
                        candidate_pos.append((pos, self.reference[pos:pos+len_frag * 3]))
                # 2. the second fragment matches the reference genome
                elif i == 1:
                    for pos in ref_pos_map[frag]:
                        candidate_pos.append((pos - len_frag, self.reference[pos - len_frag:pos+len_frag * 2]))
                # 3. the third fragment matches the reference genome
                elif i == 2:
                    for pos in ref_pos_map[frag]:
                        candidate_pos.append((pos - len_frag * 2, self.reference[pos - len_frag * 2:pos+len_frag]))
            
        return candidate_pos
    
    def find_best_pos(self, read, pos_map):
        if len(read) < 45:
            return -1, None
        
        read_length = 45
        # find the best position of the read in the reference genome
        # return the best position and the corresponding fragment from the reference genome
        candidate_pos = self.find_possible_pos(pos_map, read)

        # remove duplicates
        candidate_pos = list(set(candidate_pos))

        best_pos = -1
        best_score = float('inf')
        possible_indel_sub_loc = []
        for pos, frag in candidate_pos:
            score, possible_indel_sub_loc = self.needleman_wunsch(read[0:read_length], frag)
            if score < best_score:
                best_score = score
                best_pos = pos

        return best_pos, possible_indel_sub_loc

def main():
    reference = open('project1b-s_reference_genome.fasta', 'r')
    reference = ''.join(reference.readlines()[1:]).replace('\n', '')

    reads = open('project1b-s_with_error_paired_reads.fasta', 'r')
    reads = [line.replace('\n', '') for line in reads.readlines() if line[0] != '>']

    # read_map
    mapping = Read_mapper(reference, reads, 2, 'dp')
    mapping.read_map()



if __name__ == '__main__':
    main()
