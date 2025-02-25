from tqdm import tqdm
from collections import deque

class GenomeAssembler:
    def __init__(self, reads, headers):
        self.reads = reads
        self.headers = headers
        self.kmer_freq = self._filter_kmers(self._create_kmer_frequency_dict())
        self.graph = self._build_graph()


    def _create_kmer_frequency_dict(self):
        kmer_size = 15
        kmer_freq = {}
        for read in tqdm(self.reads, desc="Creating kmer frequency dict"):
            for i in range(len(read) - kmer_size + 1):
                kmer = read[i:i + kmer_size]
                kmer_freq[kmer] = kmer_freq.get(kmer, 0) + 1
        return kmer_freq
    
    def _filter_kmers(self, kmer_freq):
        filtered_kmer_freq = {}
        for kmer, freq in tqdm(kmer_freq.items(), desc="Filtering kmers", total=len(kmer_freq)):
            if freq >= 3:
                filtered_kmer_freq[kmer] = freq
        return filtered_kmer_freq
    
    def _build_graph(self):
        matrix = {}
        for kmer in tqdm(self.kmer_freq, desc="Building de Bruijn graph", total=len(self.kmer_freq)):
            prefix = kmer[:-1]
            suffix = kmer[1:]
            if prefix in matrix:
                matrix[prefix].append(suffix)
            else:
                matrix[prefix] = [suffix]
        return matrix


    def _find_max_non_branching_path(self):
        # calculagte the in and out degree of each node
        deg = {}
        for out_node, in_nodes in self.graph.items():
            if out_node not in deg:
                deg[out_node] = [0, len(in_nodes)]
            else:
                deg[out_node][1] += len(in_nodes)

            for in_node in in_nodes:
                if in_node not in deg:
                    deg[in_node] = [1, 0]
                else:
                    deg[in_node][0] += 1
        paths = []
        visited = set()
        # find all maximal non-branching paths, and start from unbalanced nodes
        for v in self.graph:
            if deg[v] != [1, 1]:
                if deg[v][1] > 0:
                    for w in self.graph[v]:
                        path = [v, w]
                        while deg[w] == [1, 1]:
                            if w in visited:
                                break
                            visited.add(w)
                            w = self.graph[w][0]
                            path.append(w)
                        paths.append(path)

        # add isolated cycles to the paths TODO: check if this is correct
        for v in self.graph:
            if deg[v] == [1, 1] and v not in visited:
                cycle = [v]
                visited.add(v)
                w = self.graph[v][0]
                while w != v:
                    visited.add(w)
                    cycle.append(w)
                    w = self.graph[w][0]
                paths.append(cycle)

        return paths

    def _build_contigs(self, paths):
        contigs = []
        for path in paths:
            contig = path[0]
            for kmer in path[1:]:
                contig += kmer[-1]
            contigs.append(contig)
        return contigs
    
    def _reconstruct_genome(self):
        paths = self._find_max_non_branching_path()
        contigs = self._build_contigs(paths)
        genome = ""
        for contig in contigs:
            genome += contig
        return genome

    def _map_reads_to_genome(self):
        genome = self._reconstruct_genome()
        kmer_size = 15  # same kmer size as used earlier
        genome_index = {}
        for i in range(len(genome) - kmer_size + 1):
            kmer = genome[i:i+kmer_size]
            genome_index.setdefault(kmer, []).append(i)
        
        read_positions = []
        for header, read in tqdm(zip(self.headers, self.reads), total=len(self.reads), desc="Mapping reads"):
            candidate_votes = {}
            for i in range(len(read) - kmer_size + 1):
                kmer = read[i:i+kmer_size]
                if kmer in genome_index:
                    for pos in genome_index[kmer]:
                        candidate_start = pos - i
                        if candidate_start < 0 or candidate_start > len(genome) - len(read):
                            # wrong i think
                            continue
                        candidate_votes[candidate_start] = candidate_votes.get(candidate_start, 0) + 1
            
            if candidate_votes:
                best_candidate = max(candidate_votes, key=candidate_votes.get)
            else:
                best_candidate = len(genome)
            
            read_positions.append((best_candidate, header))
        
        read_positions.sort(key=lambda x: x[0])
        return [element[1] for element in read_positions]




    def write_output(self, output_file):
        sorted_headers = self._map_reads_to_genome()
        with open(output_file, "w") as f:
            for header in tqdm(sorted_headers):
                f.write(f">{header}\n")

def load_reads(file_path):
    reads, headers = [], []
    with open(file_path, "r") as f:
        for line in tqdm(f, desc="Loading reads"):
            line = line.strip()
            if line.startswith(">"):
                headers.append(line[1:])
            else:
                reads.append(line)
    return reads, headers

def main():
    reads, headers = load_reads("project2b_reads.fasta")
    assembler = GenomeAssembler(reads, headers)
    assembler.write_output("predictions.txt")

if __name__ == "__main__":
    main()
