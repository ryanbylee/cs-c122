class Read_map():
    def __init__(self, reference, read, threshold):
        self.reference = reference
        self.read = read
        self.threshold = threshold

    def read_map(self):
        # sliding window
        # create predcited_mutations.txt

        mutation_list = open('predicted_mutations.txt', 'w')
        dup_list = set()
        read_length = 50
        for k, read in enumerate(self.read):
            if read[0] == '>':
                continue
            else:
                read = read.replace('\n', '')
                for i in range(len(self.reference) - read_length + 1):
                    window = self.reference[i:i+read_length]
                    count = 0
                    for j in range(read_length):
                        if read[j] != window[j]:
                            count += 1

                    if count <= self.threshold:
                        for j in range(read_length):

                            if read[j] != window[j]:
                                if i+j not in dup_list:
                                    dup_list.add(i + j)
                                    mutation_list.write('>S' + str(i+j) + ' ' + window[j] + ' ' + read[j] + '\n')
        mutation_list.close()



def main():
    reference = open('sample_reference_genome.fasta', 'r')

    # delete the first line

    # combine all lines into one string
    reference = reference.readlines()[1:]

    # combine all elements of the list into one string
    reference = ''.join(reference)

    # remove the newline character
    reference = reference.replace('\n', '')

    sample_no_error_single = open('sample_no_error_single_reads.fasta', 'r')

    # read_map
    mapping = Read_map(reference, sample_no_error_single, 3)
    mapping.read_map()


    # sanity check
    pred = open('predicted_mutations.txt', 'r')
    ans = open('sample_mutations.txt', 'r')

    # look for inclusion between the two files
    pred = pred.readlines()
    ans = ans.readlines()

    # print(pred)
    ans_filtered = [line for line in ans if line[0:2] == '>S']

    # print(ans_filtered)

    for line in ans_filtered:
        if line not in pred:
            print("not matched", line)


if __name__ == '__main__':
    main()
