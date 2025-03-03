import torch.nn as nn
import torch.nn.functional as F
from tqdm import tqdm
import torch



class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Conv2d(1,16, kernel_size=(21, 4), stride=1)
        self.conv2 = nn.Conv2d(16, 32, kernel_size=(5, 1), stride=1)
        self.pool = nn.MaxPool2d(kernel_size=(2, 1))
        self.fc1 = nn.Linear(32*43, 128)
        self.fc2 = nn.Linear(128, 1)

        # dimension calculation:
        '''
        initiai: (N, 1, 201, 4)
        after conv1: (N, 8, 191, 4)
        after pool: (N, 8, 95, 4)
        after conv2: (N, 16, 91, 1)
        after pool: (N, 16, 45, 1)
        after flatten: (N, 16*45)
        '''

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = self.pool(x)
        x = F.relu(self.conv2(x))
        x = self.pool(x)
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x
    


def train_valid(model, train_loader, valid_loader, criterion, optimizer, num_epochs):
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

def test(model, test_loader):
    total = 0
    correct = 0
    for inputs, labels in test_loader:
        outputs = model(inputs)
        predicted = (torch.sigmoid(outputs) > 0.5).float()
        total += labels.size(0)
        correct += (predicted.squeeze() == labels).sum().item()
    print(f'Test Accuracy: \t{100 * correct / total:.2f}%')