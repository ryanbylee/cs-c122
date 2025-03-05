from DNA_Dataset import *
from model import *
from utils import *
from torch.utils.data import DataLoader
import os
import numpy as np
import tensorboard


def main():
    torch.manual_seed(42)

    # check if GPU is available

    device_config = None
    if torch.cuda.is_available():
        device_config = 'cuda'
    elif torch.backends.mps.is_available():
        device_config = 'mps'
    else:
        device_config = 'cpu'
    device = torch.device('cpu')

    print(f'Using device: {device}')
    model = Model().to(device)
    print(model._parameters)

    weights_exist = os.path.isfile('weights.pth')
    if weights_exist:
        print('Loading model weights...')
        model.load_state_dict(torch.load('weights.pth'))
    else:
        print('Training new model...')
        dataset = DNA_Dataset("bound.fasta", "notbound.fasta", 'cpu')
        print(f'shape: {dataset.data.shape}')

        trainset, validset = torch.utils.data.random_split(dataset, [0.95, 0.05]) # 90 10 split for best res
        print(f"trainset: {len(trainset)}, validset: {len(validset)}")
        train_loader = DataLoader(trainset, batch_size=16, shuffle=True)
        valid_loader = DataLoader(validset, batch_size=16, shuffle=True)
        # test_loader = DataLoader(testset, batch_size=32, shuffle=False)
        criterion = nn.BCEWithLogitsLoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

        train_losses, valid_losses = train_valid(model, train_loader, valid_loader, criterion, optimizer, num_epochs=6)


        # test(model, test_loader)

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

    plot_losses(train_losses, valid_losses)

if __name__ == "__main__":
    main()