import torch.nn as nn
import torch.nn.functional as F
from tqdm import tqdm
import torch



class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()

        self.pool_weight = nn.Parameter(torch.tensor(0.5))

        self.seq1 = nn.Sequential(
            nn.Conv2d(1, 16, kernel_size=(21, 4), stride=1),
            nn.ReLU(),
        )
        self.seq2 = nn.Sequential(
            nn.Conv2d(16, 32, kernel_size=(5, 1), stride=1),
            nn.ReLU(),
        )
        
        self.fc1 = nn.Linear(32*43, 128)
        self.fc2 = nn.Linear(128, 1)
        # self.fc1 = nn.Linear(16*90, 128)
        # self.fc2 = nn.Linear(128, 1)

        # dimension calculation:
        '''
        for 2 conv:
        initiai: (N, 1, 201, 4)
        after conv1: (N, 16, 177, 1)
        after pool: (N, 16, 88, 1)
        after conv2: (N, 32, 84, 1)
        after pool: (N, 32, 42, 1)
        after flatten: (N, 32*42)
        '''
        
        '''
        for 1 conv:
        initial: (N, 1, 201, 4)
        after conv1: (N, 16, 181, 1)
        after pool: (N, 16, 90, 1)
        after flatten: (N, 16*90)
        '''

        '''
        for 3 conv:
        initial: (N, 1, 201, 4)
        after conv1: (N, 16, 181, 1)
        after pool: (N, 16, 90, 1)
        after conv2: (N, 32, 81, 1)
        after pool: (N, 32, 40, 1)
        after conv3: (N, 64, 36, 1)
        after pool: (N, 64, 18, 1)
        after flatten: (N, 64*18)

        '''
    def pooling(self, x):
        x = (self.pool_weight * nn.MaxPool2d(kernel_size=(2, 1), stride=2)(x)) + ((1-self.pool_weight) * nn.AvgPool2d(kernel_size=(2, 1), stride=2)(x))
        return x
    
    def forward(self, x):
        x = self.seq1(x)
        x = self.pooling(x)
        x = self.seq2(x)
        x = self.pooling(x)
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x
    


def train_valid(model, train_loader, valid_loader, criterion, optimizer, num_epochs):
    # create plot
    train_losses = []
    valid_losses = []

    # TODO: 5-fold cross validation?
    # for fold in range(5):
    #     train_loader, valid_loader = get_fold(fold)
    #     train_loader = DataLoader(trainset, batch_size=16, shuffle=True)
    #     valid_loader = DataLoader(validset, batch_size=16, shuffle=True)
    #     test_loader = DataLoader(testset, batch_size=32, shuffle=False)
    
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for _, (inputs, labels) in tqdm(enumerate(train_loader), total=len(train_loader), desc=f"Epoch {epoch + 1}/{num_epochs}"):
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels.float())
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        train_losses.append(running_loss / len(train_loader))

        print(f'Train Loss: \t\t{running_loss / len(train_loader):.4f}')

        model.eval()
        with torch.no_grad():
            valid_loss = 0.0
            correct = 0
            total = 0
            for inputs, labels in valid_loader:
                outputs = model(inputs)
                loss = criterion(outputs.squeeze(), labels.float())
                valid_loss += loss.item()
                predicted = (torch.sigmoid(outputs) > 0.5).float()
                total += labels.size(0)
                correct += (predicted.squeeze() == labels).sum().item()

            print(f'Validation Loss: \t{valid_loss / len(valid_loader):.4f}, Accuracy: {100 * correct / total:.2f}%')
            valid_losses.append(valid_loss / len(valid_loader))
    return train_losses, valid_losses

def test(model, test_loader):
    total = 0
    correct = 0
    for inputs, labels in test_loader:
        outputs = model(inputs)
        predicted = (torch.sigmoid(outputs) > 0.5).float()
        total += labels.size(0)
        correct += (predicted.squeeze() == labels).sum().item()
    print(f'Test Accuracy: \t{100 * correct / total:.2f}%')