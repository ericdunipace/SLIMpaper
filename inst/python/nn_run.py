# Import packages
import torch
from inst.python.NN import neural, train
#import torch.utils.data as data_utils

def nn_fun(D, nLayer, nLayer2):

    neural_jit = torch.jit.script(neural(D, nLayer, nLayer2))
    
    return(neural_jit)

