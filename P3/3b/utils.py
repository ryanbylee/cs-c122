import torch
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import numpy as np

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


def plot_losses(train_losses, valid_losses):
    plt.plot(train_losses, label='train')
    plt.plot(valid_losses, label='valid')
    plt.legend()
    plt.xticks(range(1, len(train_losses) + 1))
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.ylim(0, 0.4)
    plt.title('Training and Validation Loss')
    plt.show()

def plot_auc(valid_loader, model):
    pred = []
    labels_pred = []
    with torch.no_grad():
        for inputs, labels in valid_loader:
            outputs = model(inputs)
            pred.append(torch.sigmoid(outputs).cpu().numpy())
            labels_pred.append(labels.cpu().numpy())
    pred = np.concatenate(pred)
    labels_pred = np.concatenate(labels_pred)


    fpr, tpr, _ = roc_curve(labels_pred, pred)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.show()
