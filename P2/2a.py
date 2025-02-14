from tqdm import tqdm

class graph():
    def __init__(self, spectrum):
        self.graph = self.build_graph()
        self.reads = spectrum


    def build_graph(self):
        graph = {}
        for i, read in enumerate(self.reads):
            prefix = read[:-1]
            suffix = read[1:]
            if prefix not in graph:
                graph[prefix] = []
            graph[prefix].append((suffix, i))
        return graph
    
    def eulerian_path(self):
        graph = self.graph
        stack = []
        path = []
        current_node = list(graph.keys())[0]
        while stack or current_node in graph:
            if current_node in graph:
                stack.append(current_node)
                next_node = graph[current_node].pop()
                current_node = next_node
            else:
                path.append(current_node)
                current_node = stack.pop()
        path.append(current_node)
        return path[::-1]
    
    def path_to_gene(self, path):
        if not path:
            return ""
        genome = path[0]
        for i in range(1, len(path)):
            genome += path[i][-1]
        return genome
    
    def reconstruct_genome(self):
        path = self.eulerian_path()
        return self.path_to_gene(path)


if __name__ == "__main__":
    reads = open('project2_sample1_spectrum.fasta', 'r')
    reads = [line.replace('\n', '') for line in tqdm(reads.readlines(), desc='loading reads...') if line[0] != '>']

    de_bruijn_graph = graph(reads)
    output = open('predictions.txt', 'w')

    path = de_bruijn_graph.eulerian_path()

    for i in path:


    genome = de_bruijn_graph.reconstruct_genome()


