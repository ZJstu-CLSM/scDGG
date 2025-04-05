import os
import pickle
import random
import time
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import StratifiedShuffleSplit, train_test_split
from torchsummary import summary
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score
from sklearn.metrics import accuracy_score
from scipy.optimize import linear_sum_assignment as linear_assignment
from config import seed, batch_size, lr, k
from net import CAE, CAE_2, weight_init


def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True



seed_torch(seed)


class my_dataset(Dataset):
    def __init__(self, x, y):
        self.x, self.y = x, y

    def __len__(self):
        return self.y.shape[0]

    def __getitem__(self, item):
        return self.x[item], self.y[item]


class Logging(object):
    def __init__(self, name='log_PBMC-Zheng4k-ablation.txt',root = 'result'):
        self.name = os.path.join(root,name)
        self.write_format = "{time_str} " \
                            "[Epoch] {e} " \
                            "[train_loss] {train_loss:.4f} " \
                            "[test_loss] {test_loss:.4f} " \
                            "[ARI] {ARI:.4f} " \
                            "[NMI] {NMI:.4f} " \
                            "[ACC] {ACC:.4f}\n"

        with open(self.name, 'w') as f:
            f.write(f"[Seed] {seed} [Lr] {lr} [Batch] {batch_size}[k] {k}\n")

    def write_(self, *args):
        e, train_loss, test_loss, ARI, NMI, ACC = args
        data = self.write_format.format(time_str=time.strftime("%Y_%m_%d_%H:%M:%S"),
                                        e=e, train_loss=train_loss, test_loss=test_loss,
                                        ARI=ARI, NMI=NMI, ACC=ACC)
        with open(self.name, "a") as f:
            f.write(data)


def load_data():
    X, y = pickle.load(open('./data/PBMC-Zheng4k/PBMC-Zheng4k-Ablation.pkl', 'rb'))
    # X.shape: n * 1 * 40 * 40
    # y.shape: n * 1
    print(y)
    print(X.shape, y.shape)
    print(np.unique(y))
    train_x, test_x, train_y, test_y = train_test_split(X, y, test_size=0.1, random_state=42)
    train_dataset = my_dataset(train_x, train_y)
    test_dataset = my_dataset(test_x, test_y)
    all_dataset = my_dataset(X, y)

    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, drop_last=True)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, drop_last=True)
    all_dataloader = DataLoader(all_dataset, batch_size=batch_size, drop_last=False)
    return all_dataloader, train_dataloader, test_dataloader


def build_model():
    
    model = CAE_2(k=k)
    model.apply(weight_init)
    model.to()
    # print(summary(model, input_size=(1, 40, 40), device='cpu'))
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    # torch.optim.lr_scheduler.StepLR(optimizer)
    return model, optimizer


def eval_one_epoch(model: nn.Module, dataloader: DataLoader):
    model.eval()
    total_loss=0
    with torch.no_grad():
        for x, y in dataloader:
            # x,y = x.to(device),y.to(device)  //将数据移到指定的设备上
            loss = model(x)
            total_loss += loss.cpu().item()
    return total_loss / len(dataloader)


from tqdm import tqdm


def train_one_epoch(model: nn.Module, optimizer: torch.optim.Optimizer, dataloader: DataLoader):
    model.train()
    total_loss = 0
    for x, y in dataloader:
        # x,y = x.to(device),y.to(device)
        loss = model(x)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        total_loss += loss.cpu().item()

    return total_loss / len(dataloader)


def clustering(model: nn.Module, dataloader: DataLoader):
    model.eval()
    flag = 0
    with torch.no_grad():
        for x, label in tqdm(dataloader):
            # x,label=x.to(device),label.to(device)
            out = model.predict(x).cpu().numpy()
            if flag == 0:
                en_out = out
                y_true = label.numpy()
                flag = 1
            else:
                en_out = np.concatenate([en_out, out], axis=0)
                y_true = np.concatenate([y_true, label], axis=0)
    # centroids, inertia, y_pred = Kmeans(pre, y, k, seed=seed)
    nmi, ari, acc = Kmeans(en_out, y_true, k, seed=seed)
    return nmi, ari, acc


def Kmeans(pre_res, y, n_clusters, n_init=40, max_iter=300, seed=2022):
    kmeans_model = KMeans(init='k-means++', n_clusters=n_clusters, n_init=n_init, max_iter=max_iter, random_state=seed)
    y_pred = kmeans_model.fit_predict(pre_res)
    D = 1.0 / euclidean_distances(pre_res, kmeans_model.cluster_centers_, squared=True)
    D **= 2.0 / (2 - 1)
    D /= np.sum(D, axis=1)[:, np.newaxis]
    centroids = kmeans_model.cluster_centers_.T

    nmi = normalized_mutual_info_score(y, y_pred)
    ari = adjusted_rand_score(y, y_pred)
    acc = cal_acc(y, y_pred)
    print(f'nmi: {nmi:.4f} ari: {ari:.4f} acc:{acc:.4f}')
    
    y_pred = pd.DataFrame(y_pred)
    y_pred.to_csv("./data/PBMC-Zheng4k/cluster_PBMC-Zheng4k-Ablation.csv")

    # return centroids, kmeans_model.inertia_, y_pred
    return nmi, ari, acc

from scipy.optimize import linear_sum_assignment as linear_assignment_


# from sklearn.utils import linear_assignment_
def cal_acc(y_true, y_pred):
    """
    Calculate clustering accuracy. Require scikit-learn installed
    # Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`
    # Return
        accuracy, in [0,1]
    """
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    ind = linear_assignment(w.max() - w)
    return sum([w[i, j] for i, j in zip(ind[0], ind[1])]) * 1.0 / y_pred.size
