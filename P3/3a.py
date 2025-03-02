from tqdm import tqdm
import random

class PeakSummitPredictor:
    def __init__(self, bounds_path, seqs_path, motif_len=21, max_iter=100, num_starts=10):
        self.seed = 42
        self.bounds = self._load_fasta(bounds_path)
        self.seqs = self._load_fasta(seqs_path)
        self.max_iter = max_iter
        self.num_inits = num_starts
        self.motif_length = motif_len

        self.pwm = self._discover_pwm()



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
    

    def _discover_pwm(self):
        nuc_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
        pwms = []
        for i in range(self.num_inits):
            tqdm.write(f"Starting random initialization {i + 1}/{self.num_inits}")
            positions = [random.randint(0, len(self.bounds[0]) - self.motif_length) for _ in range(len(self.bounds))]
            pwm = self._create_pwm(positions, nuc_to_index)

            for j in tqdm(range(self.max_iter), desc="optimizing PWM"):
                updated_pos = []
                for bound in self.bounds: 
                    best_pos, best_score = self._find_best_position(bound, pwm, nuc_to_index)
                    updated_pos.append(best_pos)
                if updated_pos == positions:
                    tqdm.write(f"Converged after {j} iterations")
                    break
                positions = updated_pos
                pwm = self._create_pwm(positions, nuc_to_index)
            if j == self.max_iter - 1:
                print('pwm not converged')
            pwm_score = self._calculate_score(pwm, nuc_to_index)
            pwms.append((pwm, pwm_score))
        best_pwm, _ = max(pwms, key=lambda x: x[1])
        return best_pwm

    
    def _find_best_position(self, seq, pwm, nuc_to_index):
        best_pos = -1
        best_score = float("-inf")
        # i is the start position of the subsequence
        for i in range(len(seq) - self.motif_length + 1):
            subseq = seq[i : i + self.motif_length]
            score = 0
            for j, nuc in enumerate(subseq):
                score += pwm[nuc_to_index[nuc]][j]
            if score > best_score:
                best_score = score
                best_pos = i
        return best_pos, best_score
    
    def _calculate_score(self, pwm, nuc_to_index):
        # calculate average likelihood
        total_score = 0
        for bound in self.bounds:
            _, best_score = self._find_best_position(bound, pwm, nuc_to_index)
            total_score += best_score
        return total_score / len(self.seqs)
    
    def _create_pwm(self,positions, nuc_to_index, motif_length=21, ):
        counts = [[0] * motif_length for _ in range(4)]
        
        for seq, start in zip(self.bounds, positions):
            subseq = seq[start : start + motif_length]
            for i, nuc in enumerate(subseq):
                counts[nuc_to_index[nuc]][i] += 1

        counts = [[count + 1 for count in row] for row in counts]

        pwm = [[count / len(positions) for count in row] for row in counts]

        return pwm

    def predict_peak_summits(self):
        nuc_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}

        predictions = []
        for i, seq in tqdm(enumerate(self.seqs), desc="Predicting peak summits", total=len(self.seqs)):
            best_pos, _ = self._find_best_position(seq, self.pwm, nuc_to_index)

            if best_pos < 31:
                best_pos = 31
            elif best_pos > 171:
                best_pos = 171
            predictions.append((i + 1, best_pos + 1))

        # write output
        with open("predictions.txt", "w") as f:
            for i, pos in predictions:
                f.write(f"seq{i},{pos}\n")




def main():
    predictor = PeakSummitPredictor("boundcentered.fasta", "boundrandomoffset.fasta")
    predictor.predict_peak_summits()


if __name__ == "__main__":
    main()



      