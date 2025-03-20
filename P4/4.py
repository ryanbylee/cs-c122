import sys

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.read().strip().splitlines()
    obs_chars = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(obs_chars)

def forward(obs, states, trans, emiss):
    T = len(obs)
    alpha = [{} for _ in range(T)]
    scale = [0.0 for _ in range(T)]
    init_prob = 1.0 / len(states)
    for s in states:
        alpha[0][s] = init_prob * emiss[s][obs[0]]
    scale[0] = sum(alpha[0][s] for s in states)
    if scale[0] == 0:
        sys.exit("Error: scaling factor at time 0 is zero.")
    for s in states:
        alpha[0][s] /= scale[0]
    for t in range(1, T):
        for s in states:
            alpha[t][s] = sum(alpha[t-1][prev] * trans[prev][s] for prev in states) * emiss[s][obs[t]]
        scale[t] = sum(alpha[t][s] for s in states)
        if scale[t] == 0:
            sys.exit(f"Error: scaling factor at time {t} is zero.")
        for s in states:
            alpha[t][s] /= scale[t]
    return alpha, scale

def backward(obs, states, trans, emiss, scale):
    T = len(obs)
    beta = [{} for _ in range(T)]
    for s in states:
        beta[T-1][s] = 1.0
    for t in range(T - 2, -1, -1):
        for s in states:
            beta[t][s] = sum(trans[s][s_next] * emiss[s_next][obs[t+1]] * beta[t+1][s_next] for s_next in states)
            beta[t][s] /= scale[t+1]
    return beta

def baum_welch(obs, states, alphabet, trans, emiss, iterations):
    T = len(obs)
    for it in range(iterations):
        alpha, scale = forward(obs, states, trans, emiss)
        beta = backward(obs, states, trans, emiss, scale)
        gamma = [ { s: 0.0 for s in states } for _ in range(T) ]
        xi = [ { (i, j): 0.0 for i in states for j in states } for _ in range(T - 1) ]
        for t in range(T - 1):
            denom = 0.0
            for i in states:
                for j in states:
                    denom += alpha[t][i] * trans[i][j] * emiss[j][obs[t+1]] * beta[t+1][j]
            if denom == 0:
                sys.exit(f"Error: Denominator zero at time {t} when computing xi.")
            for i in states:
                gamma[t][i] = 0.0
                for j in states:
                    xi_val = (alpha[t][i] * trans[i][j] * emiss[j][obs[t+1]] * beta[t+1][j]) / denom
                    xi[t][(i, j)] = xi_val
                    gamma[t][i] += xi_val
        t = T - 1
        denom = sum(alpha[t][s] * beta[t][s] for s in states)
        for s in states:
            gamma[t][s] = (alpha[t][s] * beta[t][s]) / denom if denom != 0 else 0.0
        new_trans = { i: { j: 0.0 for j in states } for i in states }
        for i in states:
            numer = sum(gamma[t][i] for t in range(T - 1))
            for j in states:
                numer_ij = sum(xi[t][(i, j)] for t in range(T - 1))
                new_trans[i][j] = numer_ij / numer if numer != 0 else 0.0
        new_emiss = { i: { symbol: 0.0 for symbol in alphabet } for i in states }
        for i in states:
            numer = sum(gamma[t][i] for t in range(T))
            for symbol in alphabet:
                numer_sym = sum(gamma[t][i] for t in range(T) if obs[t] == symbol)
                new_emiss[i][symbol] = numer_sym / numer if numer != 0 else 0.0
        trans, emiss = new_trans, new_emiss
    return trans, emiss

def compute_gamma(obs, states, trans, emiss):
    alpha, scale = forward(obs, states, trans, emiss)
    beta = backward(obs, states, trans, emiss, scale)
    T = len(obs)
    gamma = []
    for t in range(T):
        gamma_t = {}
        norm = sum(alpha[t][s] * beta[t][s] for s in states)
        for s in states:
            gamma_t[s] = (alpha[t][s] * beta[t][s]) / norm if norm != 0 else 0.0
        gamma.append(gamma_t)
    return gamma

def main():
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "input.fasta"
    obs = read_fasta(input_file)
    alphabet = ['n', 'x', 'y', 'z']
    states = ['G', 'N']
    trans = {
        "G": {"G": 0.9, "N": 0.1},
        "N": {"G": 0.1, "N": 0.9}
    }
    emiss = {
        "G": {"n": 0.1, "x": 0.3, "y": 0.3, "z": 0.3},
        "N": {"n": 0.7, "x": 0.1, "y": 0.1, "z": 0.1}
    }
    iterations = 10
    trans, emiss = baum_welch(obs, states, alphabet, trans, emiss, iterations)
    gamma = compute_gamma(obs, states, trans, emiss)
    gene_scores = { s: emiss[s]['x'] + emiss[s]['y'] + emiss[s]['z'] for s in states }
    gene_state = max(gene_scores, key=gene_scores.get)
    gene_posteriors = [ gamma[t][gene_state] for t in range(len(obs)) ]
    top_count = 50000
    sorted_indices = sorted(range(len(gene_posteriors)), key=lambda i: gene_posteriors[i])
    top_indices = sorted_indices[-top_count:]
    top_indices_sorted = sorted(top_indices)
    output = open("output.txt", "w")
    for idx in top_indices_sorted:
        output.write(f"{idx + 1}\n")
if __name__ == "__main__":
    main()
