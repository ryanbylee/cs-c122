import torch

def load_fasta(file_path):
        print(f"loading {file_path}...")
        seqs = []
        f = open(file_path, 'r').readlines()
        i = 0 
        while i < len(f):
            line = f[i].strip()
            if line.startswith(">") or line == "":
                i += 1
                continue
            else:
                seq_parts = ""
                for j in range(1, 5):
                    seq_parts += f[i+j].strip()
                seqs.append(line + seq_parts)
                i+=5
        return seqs

def one_hot_encode(seq):
    nuc_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
    one_hot_seq = torch.zeros((len(seq), 4))
    for j, nuc in enumerate(seq):
        one_hot_seq[j, nuc_to_index[nuc]] = 1
    return one_hot_seq

