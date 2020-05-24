# Import packages
import torch

def jit_fun(D, nLayer, nLayer2):
    from NN.model import neural
    
    neural_jit = torch.jit.script(neural(D, nLayer, nLayer2))
    
    return(neural_jit)

