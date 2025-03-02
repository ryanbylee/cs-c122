from tqdm import tqdm
from collections import deque

class DeBruijnGraph:
    def __init__(self, spectrum, headers):
        self.reads = spectrum
        self.headers = headers
        self.graph = self._build_graph()
        self.path = self._find_eulerian_path()

    def _build_graph(self):
        g = {}
        for i, kmer in tqdm(enumerate(self.reads), total=len(self.reads)):
            prefix = kmer[:-1]
            suffix = kmer[1:]
            if prefix not in g:
                g[prefix] = []
            if suffix not in g:
                g[suffix] = []
            g[prefix].append(suffix)
        return g

    def _find_eulerian_path(self):
        in_deg, out_deg = {}, {}
        nodes = set(self.graph.keys())
        for node, neighbors in tqdm(self.graph.items(), desc="in/out deg"):
            out_deg[node] = out_deg.get(node, 0) + len(neighbors)
            for neighbor in neighbors:
                in_deg[neighbor] = in_deg.get(neighbor, 0) + 1
                nodes.add(neighbor)
                if neighbor not in self.graph:
                    self.graph[neighbor] = []
        for node in nodes:
            if node not in out_deg:
                out_deg[node] = 0
            if node not in in_deg:
                in_deg[node] = 0
        start_node = next(iter(self.graph))
        for node in nodes:
            if out_deg.get(node, 0) - in_deg.get(node, 0) == 1:
                start_node = node
                break
        stack = [start_node]
        path = deque()
        while stack:
            node = stack[-1]
            if node not in self.graph:
                self.graph[node] = []
            if self.graph[node]:
                stack.append(self.graph[node].pop(0))
            else:
                path.appendleft(stack.pop())
        return list(path)

    def reconstruct_genome(self):
        genome = self.path[0]
        for kmer in tqdm(self.path[1:]):
            genome += kmer[-1]
        return genome

    def map_reads_to_genome(self):
        genome = self.reconstruct_genome()
        ordered_headers = []
        for header, read in tqdm(zip(self.headers, self.reads), total=len(self.reads)):
            if read in genome:
                ordered_headers.append((genome.index(read), header))
        ordered_headers.sort()
        return [h[1] for h in ordered_headers]

    def write_output(self, output_file):
        sorted_headers = self.map_reads_to_genome()
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
    reads, headers = load_reads("project2a_spectrum.fasta")
    db_graph = DeBruijnGraph(reads, headers)
    db_graph.write_output("predictions.txt")

if __name__ == "__main__":
    main()
