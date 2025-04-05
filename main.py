import numpy as np
import torch
from torchsummary import summary
import matplotlib.pyplot as plt



from net import CAE, weight_init
from utils import load_data,build_model,train_one_epoch,eval_one_epoch,clustering,Logging
if __name__ == '__main__':
    device = torch.device("cuda" if  torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    all_dataloader,train_dataloader, test_dataloader = load_data()

    model,optimizer = build_model()

    Log = Logging(name='PBMC-Zheng4k-Ablation.txt')
    best_model = 0
    train_losses = []  # 存储每个 epoch 的训练损失
    test_losses = []   # 存储每个 epoch 的测试损失
    acc_fin = 0
    for e in range(200):
        train_loss = train_one_epoch(model, optimizer, train_dataloader)
        test_loss = eval_one_epoch(model,test_dataloader)

        # 将损失值添加到列表中
        train_losses.append(train_loss)
        test_losses.append(test_loss)
        
        print(f'Epoch: {e+1} train: {train_loss:.4f} test: {test_loss:.4f}')
        nmi,ari,acc = clustering(model,all_dataloader)
        print('-'*50)
        Log.write_(e+1,train_loss,test_loss,ari,nmi,acc)
        if acc > best_model:
            acc_fin =acc
            print(f"best:{acc}")
            best_model = acc
            torch.save(model.state_dict(),'checkpoint/PBMC-Zheng4k-Ablation.pth')
    print(f"acc_fin:{acc_fin}")       
    # 绘制训练和测试损失的变化曲线
    plt.figure(figsize=(10, 6))
    plt.plot(train_losses, label='Train Loss')
    plt.plot(test_losses, label='Test Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Training and Testing Loss Over Epochs')
    plt.legend()
    plt.grid(True)
    plt.savefig('./data/PBMC-Zheng4k/loss_curves_PBMC-Zheng4k-Ablation.png')  # 保存图像
    plt.show()