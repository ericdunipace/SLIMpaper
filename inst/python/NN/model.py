# Import packages
import numpy as np
import scipy
import torch
import torch.autograd as autograd
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
#from itertools import compress,product, combinations,permutations, islice,chain

# Model Module
class neural(torch.nn.Module):
    def __init__(self, D: int, nLayer: int, nLayer2: int ) -> None:
        super(neural, self).__init__()
        self.D = D
        self.nLayer = nLayer
        self.nLayer2 = nLayer2 # intermediate layer sizes
        
        self.layer1 = nn.Linear(self.D, self.nLayer)
        self.layer2 = nn.Linear(self.nLayer, self.nLayer2)
        self.layer3 = nn.Linear(self.nLayer2, self.nLayer2)
        #self.layer4 = nn.Linear(self.nLayer2,self.nLayer2)
        self.layer5 = nn.Linear(self.nLayer2, 1)
        self.loss = nn.BCEWithLogitsLoss(reduction = 'mean')
    
    def forward(self, y: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        
        eta  =   self.predict(x)
        loss =   self.loss(eta, y)
        
        return(loss)
    
    @torch.jit.export
    def predict(self, x: torch.Tensor) -> torch.Tensor:
        
        eta  =   F.elu(self.layer1(x))
        #eta  =   F.dropout(eta)
        eta  =   F.elu(self.layer2(eta))
        #eta  =   F.dropout(eta)
        eta  =   F.elu(self.layer3(eta))
        #eta  =   F.dropout(eta)
        #eta  =   F.elu(self.layer4(eta))
        eta  =   F.elu(self.layer5(eta))
        
        return(eta)


# training function
def train(model: torch.nn.Module,
train_iter: torch.utils.data.dataloader.DataLoader,
test_iter,
num_epoch: int,
lr: float,
lam = 0.0,
test=False,
verbose=False):
    
    @torch.no_grad()
    def calc_loss(model, test_iter):
        Nt = 0
        N = 0
        test_ll = 0.0
        count = 0
        for test_batch in enumerate(test_iter):
            Xtest = test_batch[1][0].float()
            Ytest = test_batch[1][1].float()
            Nt = len(Ytest)
            test_ll = model.forward(Ytest, Xtest) * Nt / (N + Nt) + test_ll * N/(N + Nt)
            N += Nt
        return(test_ll.data.numpy())
    
    @torch.no_grad()
    def calc_acc(model, test_iter):
        N      = 0.0
        count  = 0
        acc    = 0.0
        eta    = []
        pseudo = []
        for test_batch in enumerate(test_iter):
            Xtest  = test_batch[1][0].float()
            Ytest  = test_batch[1][1].float()
            N     += len(Ytest)
            eta    = model.predict(Xtest)
            pseudo = torch.round(torch.sigmoid(eta))
            acc   += torch.sum(Ytest == pseudo)
        acc   /= N
        return(acc.data.numpy())
    
    # NN model set-up
    model.train() # set train to true, just in case
    
    # optimizer
    opt      = optim.Adam(model.parameters(), lr=lr, weight_decay=lam)
    
    # holders for iterations
    loss     = []
    dat_ll   = np.zeros(num_epoch)
    test_ll  = np.zeros(num_epoch)
    dat_acc  = np.zeros(num_epoch)
    test_acc = np.zeros(num_epoch)
    
    # iterate through training data
    for epoch in range(num_epoch):
        N = 0.0
        Ntest = 0.0
        count=0
        for t_batch in enumerate(train_iter):
            count += 1
            
            # get data
            Yt     = Variable(t_batch[1][1].float(), requires_grad=False)
            Xt     = Variable(t_batch[1][0].float(), requires_grad=False)
            #Nt = len(Yt)
            
            # zero gradients
            model.zero_grad()
            
            #calculate forward (likelihood)
            loss   = model.forward(Yt,Xt)
            
            #calculate gradients
            loss.backward()
            
            #update parameters
            opt.step()
            #dat_ll[epoch] = loss.data.numpy() * Nt/(N+Nt) + dat_ll[epoch] * N/(N+Nt)
            #N += Nt
            #print("Epoch Num: ", epoch+1,", Train: ", count, ", Current loss: ", loss.data.numpy(), end='\r')
        #run test
        if verbose:
            if test:
                model.eval()
                # with torch.no_grad():
                dat_ll[epoch]   = loss.data.numpy() #calc_loss(model, train_iter)
                test_ll[epoch]  = calc_loss(model, test_iter)
                dat_acc[epoch]  = torch.sum(Yt == torch.round(torch.sigmoid(model.predict(Xt))))/float(len(Yt)) #calc_acc(model, train_iter)
                test_acc[epoch] = calc_acc(model, test_iter)
                model.train()
                print("Epoch Num:", epoch+1,", Cur. Avg Loss:", round(dat_ll[epoch],3),
                ", Cur. Train acc:", round(dat_acc[epoch], 3),
                ", Test loss:", round(test_ll[epoch],3), ", Test acc:", round(test_acc[epoch],3), flush = True)
            else:
                dat_ll[epoch] = loss.data.numpy()
                print("Epoch Num:", epoch+1,", Current Avg Loss:", round(dat_ll[epoch],3), flush = True)
    #set to eval mode
    model.eval()
    return(dat_ll, test_ll, model)

