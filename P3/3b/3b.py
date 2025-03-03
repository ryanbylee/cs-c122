from DNA_Dataset import *
from model import *
from utils import *
from torch.utils.data import DataLoader
import os
import numpy as np

def main():
    torch.manual_seed(42)

    # check if GPU is available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Using device: {device}')
    model = Model().to(device)

    exists_weights = os.path.isfile('weights.pth')
    if exists_weights:
        print('Loading model weights...')
        model.load_state_dict(torch.load('weights.pth'))
    else:
        print('Training new model...')
        dataset = DNA_Dataset("bound.fasta", "notbound.fasta")
        print(f'shape: {dataset.data.shape}')

        trainset, validset, testset = torch.utils.data.random_split(dataset, [0.6, 0.2, 0.2])
        train_loader = DataLoader(trainset, batch_size=32, shuffle=True)
        valid_loader = DataLoader(validset, batch_size=32, shuffle=True)
        test_loader = DataLoader(testset, batch_size=32, shuffle=False)
        criterion = nn.BCEWithLogitsLoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

        train_valid(model, train_loader, valid_loader, criterion, optimizer, num_epochs=10)
        test(model, test_loader)

        # save weights
        torch.save(model.state_dict(), 'weights.pth')

    model.eval()


    # load test data
    test_data = load_fasta("test.fasta")
    # given test data, output the sequence numbers that has predicted label 1
    test_data = [one_hot_encode(seq) for seq in test_data]
    test_data = torch.stack(test_data, dim=0).unsqueeze(1)
    test_data = test_data.to(device)
    outputs = model(test_data)
    predicted = (torch.sigmoid(outputs) > 0.5).float()
    predicted = predicted.cpu().numpy()

    with open("predictions.txt", "w") as f:
        for seq_num in tqdm(np.where(predicted == 1)[0], desc="Writing predictions"):
            f.write(f"seq{seq_num + 1}\n")


if __name__ == "__main__":
    main()