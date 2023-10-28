import torch 
import torchdiffeq
import torch.nn as nn
import math
import numpy as np
from torchdiffeq import odeint_adjoint
from torchdiffeq import odeint

import logging

# save output in a file called test.out
logging.basicConfig(filename='param_estim.out', encoding='utf-8', level=logging.INFO)


#phi(a, gamma_a) = a ./ (a + gamma_a)

class model(nn.Module):
    def __init__(self, params, device):
        super().__init__()
        self.params = nn.Parameter(params)
        
    def phi(self, u, γu):
        return u / (u+γu)
        
    def forward(self, t, u):
        βa, βB, βE1, βE2, βh1, βh2, βh3, βM1, βM2, βp, γa, γB, γh, γp, μaE, μaM, μhM, μpB, μpE, q = self.params
        B, E, M, a, h, p = u
        
        
        nextB = βB*self.phi(p,γp)*B - q*B
        nextE = (βE1*self.phi(a, γa) + βE2*(1-self.phi(B, γB))*self.phi(p, γp))*E - q*E
        nextM = (βM1*self.phi(a, γa) + βM2*self.phi(h, γh))*M - q*M

        nexta = βa*self.phi(p, γp)*B - q*a - (μaE*E + μaM*M)*self.phi(a, γa)
        nexth = βh1*self.phi(a, γa)*E + βh2*self.phi(p, γp)*B - q*h- μhM*self.phi(h, γh)*M + βh3*(1-self.phi(B, γB))*self.phi(p, γp)*E 
            
        nextp = βp*q*(math.cos(t)+1)**3 - q*p - (μpB*B + μpE*E)*self.phi(p, γp)
        
        return torch.stack([nextB,nextE, nextM, nexta, nexth, nextp])
    
u0 = torch.tensor([0.0004706, 0.0004706, 0.0004706, 9.7079, 7.9551, 32.061])
tspan = (0.0,500.0);
p = torch.nn.Parameter(torch.tensor([1e3, # betaa
    0.35, # betab
    0.4, # betaE1
    0.2, # betaE2
    75, #betah1
    150, #betah2
    10,  #betah3                   
    1.0, #betaM1
    0.9, #betaM2
    100, #betap
    100, #gammaa
    1000, # gammaB
    5000, #gammah
    100, # gammap
    2500, # muaE
    5000, # muaM
    4000, # muhM
    2000, #mupB
    5000, #mupE
    0.1]),requires_grad = True) #0.1 # q

device = 'cpu'
diffeq = model(p, device)


adjoint = True
ode_integrator = odeint_adjoint if adjoint else odeint


n_points = 500
tmax = 400
#X = ode_integrator(diffeq, u0, torch.linspace(0, tmax, n_points, device=device))


# truth = u0; t_start = 400; t_end = 500; weight = u0
def loss(sol, truth, t_start, t_end, weight):
    """ loss computes a weighted mean squared error loss between a prediction (averaged across timepoints t_start to t_end) and ground truth weighted by a weight vector
    truth (Vector): vector containing the ground truth
    pred (Matrix): matrix with shape (variables, timepoints)
    t_start (integer): first timepoint to pull to compute an average across time
    t_end (integer): last timepoint to pull to compute an average across time
    """
    #print((truth-torch.mean(sol[t_start:t_end,:], dim=0))**2 /weight)
    #return sum((truth-torch.mean(sol[t_start:t_end,:], dim=0))**2 /weight )#(u0-mean(sol[:,t_start:t_end], dims=2))./ weight
    s = (sum(truth[:3]) - sum(torch.mean(sol[t_start:t_end,:3], dim =0)) )**2 / sum(weight[:3])
    return s + sum((truth[3:]-torch.mean(sol[t_start:t_end,3:], dim=0))**2 /weight[3:] ) #+ torch.var(sol[t_start:t_end,5])/weight[5]#(u0-mean(sol[:,t_start:t_end], dims=2))./ weight

n_epochs = 1000

optimizer = torch.optim.Adam([{'params':diffeq.parameters(), 'lr':1e-2}])
lambda1 = lambda epoch: (epoch-9)**(-1) if epoch >=100 else 1
scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=[lambda1])

training_loss = []

for epoch in range(n_epochs):
    
    logging.info("epoch {}".format(epoch))
    # zero out gradient information
    optimizer.zero_grad()
    
    # get ODE solution for timepoints 0 to 500, 1000 discrete points
    X = ode_integrator(diffeq, u0, torch.linspace(0, 300, 500, device=device))
    
    # compute the loss 
    l = loss(X,u0,300,500, u0)
    logging.info("current loss {}".format(l))
    
    # store loss information, do not track gradient information for this operation
    with torch.no_grad():
        training_loss.append(l.detach().cpu())
    
    # compute gradient information
    l.backward()
    
    # update parameters based on gradient information
    optimizer.step()
    
    # update scheduler
    scheduler.step()
    
    torch.save(diffeq.params,"params9_27.pt")
    logging.info("Saved current parameter settings to params9_27.pt")
