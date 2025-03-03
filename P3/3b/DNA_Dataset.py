import torch
from torch.utils.data import Dataset
from utils import *


class DNA_Dataset(Dataset):
    def __init__(self, bounds_path, notBounds_path):
        self.data, self.labels = self._prepare_data(bounds_path, notBounds_path)

    def _prepare_data(self, bounds_path, notBounds_path):
        bounds = load_fasta(bounds_path)
        notBounds = load_fasta(notBounds_path)

        bounds = [one_hot_encode(seq) for seq in bounds]
        notBounds = [one_hot_encode(seq) for seq in notBounds]

        data = torch.stack(bounds + notBounds, dim=0).unsqueeze(1)
        labels = torch.cat((torch.ones(len(bounds)), torch.zeros(len(notBounds))))


        return data, labels
    
    
    
    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]