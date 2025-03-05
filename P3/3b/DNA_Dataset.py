import torch
from torch.utils.data import Dataset
from utils import *


class DNA_Dataset(Dataset):
    def __init__(self, bounds_path, notBounds_path, device_config):
        self.device = device_config
        self.data, self.labels = self._prepare_data(bounds_path, notBounds_path)

    
    def _reverse_complement(self, seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement[base] for base in reversed(seq))

    def _prepare_data(self, bounds_path, notBounds_path):
        bounds = load_fasta(bounds_path)
        bounds_rev = [self._reverse_complement(seq) for seq in bounds]
        bounds = bounds + bounds_rev

        notBounds = load_fasta(notBounds_path)
        notBounds_rev = [self._reverse_complement(seq) for seq in notBounds]
        notBounds = notBounds + notBounds_rev

        bounds = [one_hot_encode(seq) for seq in bounds]
        notBounds = [one_hot_encode(seq) for seq in notBounds]

        data = torch.stack(bounds + notBounds, dim=0).unsqueeze(1)
        labels = torch.cat((torch.ones(len(bounds)), torch.zeros(len(notBounds))))


        return data.to(self.device), labels.to(self.device)
    
    
    
    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]