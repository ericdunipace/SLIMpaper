# Import packages
import numpy as np
import scipy
import torch
import torch.autograd as autograd
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data as data_utils
from itertools import compress,product, combinations,permutations, islice,chain


train_pd = pd.read_csv("../Data/train_melt.csv", delimiter=',')
test_pd = pd.read_csv("../Data/test_melt.csv", delimiter=',')
Std = np.genfromtxt("../Data/PM/clean/pm_stnd_param.csv", delimiter=',')
stdmu = Std[1,:]
stdsd = Std[2,:]


# In[5]:

Yt = pd.DataFrame.as_matrix(train_pd.iloc[:,0])
Ytst = pd.DataFrame.as_matrix(test_pd.iloc[:,0])


# In[6]:

train_features = torch.FloatTensor(pd.DataFrame.as_matrix(train_pd.iloc[:,1:]))
train_target = torch.FloatTensor(Yt)
test_features = torch.FloatTensor(pd.DataFrame.as_matrix(test_pd.iloc[:,1:]))
test_target = torch.FloatTensor(Ytst)


# In[7]:

train_data = data_utils.TensorDataset(train_features, train_target)
train_loader = data_utils.DataLoader(train_data, batch_size=100000, shuffle=True)
test_data = data_utils.TensorDataset(test_features, test_target)
test_loader = data_utils.DataLoader(test_data, batch_size=100000, shuffle=False)


# In[8]:

D = train_features.shape[1]-1
nLayer = D * 10
num_epoch=20
p=1
lr=0.0001
test=True


train_l, test_l, mod = train(neural, train_loader, test_loader, D, nLayer, num_epoch, lr, test)


test_mse = mod.test(test_loader, stdmu, stdsd)


print("MSE NN: ", test_mse)


# In[26]:

mselist = [test_mse]
outpath="../Output/MSE_NN.csv"
with open(outpath, "w") as file:
    writer = csv.writer(file, delimiter=',')
    writer.writerow(mselist)

mselist = [train_l.numpy(),test_l.numpy()]
outpath="../Output/NN_MSE_list.npy"
# with open(outpath, "w") as file:
#     writer = csv.writer(file, delimiter=',')
#     writer.writerow(mselist)
np.save(outpath,mselist)


# In[27]:

pred_obs = np.zeros([366,G])
predicted = []
Ytmat = np.genfromtxt('../Data/PM/clean/pm_2.5_2004-2004_final.csv', delimiter=',')[1:,:]
for test_batch in enumerate(test_loader):
    Xtest = Variable(test_batch[1][0][:,1:].float(), volatile=True)
    GroupTest = Variable(test_batch[1][0][:,0].long()-1, volatile=True)
    Ytest = Variable(test_batch[1][1].float(), volatile=True)
    predicted.append(mod.predict(Xtest, GroupTest).data.numpy())


# In[28]:

pred_out = np.vstack(predicted)
count = 0
for i in range(pred_obs.shape[0]):
    for j in range(pred_obs.shape[1]):
        pred_obs[i,j] = pred_out[count]
        count += 1


# In[29]:

pred_adjust = pred_obs.copy()
for i in range(pred_adjust.shape[0]):
    pred_adjust[i,:] = pred_adjust[i,:] *stdsd + stdmu


# In[30]:

outpath="../Output/NN_pred.npy"
np.save(outpath, pred_adjust)


# In[ ]:




