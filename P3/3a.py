from tqdm import tqdm

class PeakSummitPredictor:
    def __init__(self, bounds_path, seqs_path):
        self.bounds = self._load_fasta(bounds_path)
        self.seqs = self._load_fasta(seqs_path)
        self.pwm = self._create_pwm()


    def _load_fasta(self, file_path):
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
    
    def _create_pwm(self, motif_length=21):

        pwm = [[0 for _ in range(4)] for _ in range(motif_length)]
        nuc_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
        for seq in tqdm(self.bounds, desc="Creating PWM"):
            mid = len(seq) // 2
            for i in range(mid-5, mid+5):
                nucleotide = seq[i]
                pwm[i-mid+5][nuc_to_index[nucleotide]] += 1

        for i in range(len(pwm)):
            for j in range(len(pwm[i])):
                pwm[i][j] += 1
        
        for i in range(len(pwm)):
            total = sum(pwm[i])
            if total > 0:
                pwm[i] = [count / total for count in pwm[i]]

        return pwm




def main():
    predictor = PeakSummitPredictor("boundcentered.fasta", "boundrandomoffset.fasta")



if __name__ == "__main__":
    main()
      