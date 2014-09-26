#!/usr/bin/env python

import scipy.integrate as spi
import scipy.optimize
import numpy as np
import pylab as pl
import warnings
from matplotlib.font_manager import FontProperties as fmp
import sys
import pdb
import itertools


def make_inputs(N,nlocs,R0s,alpha=0.0,gamma=1/3.0,alpha_hs_mult=[0.0]):
    """
    returns array of inputs for main SIR model
    """
    rc =  make_inputs_helper(
        S0=N-np.array([1]*nlocs),
        I0=np.array([1]*nlocs),
        R0=np.array([0]*nlocs),
        CI0=np.array([0]*nlocs),
        V0=np.array([0]*nlocs),
        R0s=R0s,
        alpha=alpha,
        gamma=gamma,
        alpha_hs_mult=alpha_hs_mult
        )
    return rc

def make_inputs_VE(N,nlocs,R0s,alpha=0.0,theta=0.866,gamma=1/3.0,alpha_hs_mult=[0.0]):

    """
    returns array of inputs for the model that allows for an imperfect all-or-nothing vaccine
    """
    rc =  make_inputs_helper_VE(
        S0=N-np.array([1]*nlocs),
        SV0=np.array([0]*2),
        I0=np.array([1]*nlocs),
        R0=np.array([0]*nlocs),
        CI0=np.array([0]*nlocs),
        V0=np.array([0]*nlocs),
        R0s=R0s,
        alpha=alpha,
        theta=theta,
        gamma=gamma,
        alpha_hs_mult=alpha_hs_mult
        )
    return rc

def make_inputs_multiepidemic(N,nlocs,R0s,alpha=0.0,gamma=1/3.0,
                              alpha_hs_mult=[0.0]):
    """
    returns array of inputs for the multi-epidemic model (run for mulitple years
    to get a reasonable susceptibile landscape)
    """
    rc =  make_inputs_helper(
        S0=N-np.array([1]*nlocs),
        I0=np.array([1]*nlocs),
        R0=np.array([0]*nlocs),
        CI0=np.array([0]*nlocs),
        V0=np.array([0]*nlocs),
        R0s=R0s,
        alpha=alpha,
        gamma=gamma,
        alpha_hs_mult=alpha_hs_mult
        )
    return rc

def epidemic_going(run,vac_end_threshold):
    """determines if epidemic is over based on incidence threshold
    input:
    1. run = output from odeint
    2. vac_end_threshold = minimum number of new cases in a time step
    returns:
    boolean for whether epidemic is still raging if True then still going
    """
    import scipy.integrate as spi
    import numpy as np
    import pylab as pl
    import warnings
    from matplotlib.font_manager import FontProperties as fmp
    import sys
    import pdb
    import vacfunctions as vf

    # get cumulative incidence
    ci = vf.get_ci_columns(run)
    #print ci

    cases_per_step = np.array([np.diff(ci[:,x]) for x in range(ci.shape[1])])
    cases_per_step = np.array([sum(locs) for locs in cases_per_step.T])

    # figure out when the epidemic ended
    epidemic_end = ((cases_per_step < vac_end_threshold).nonzero())[0]

    return len(epidemic_end) == 0


## makes inputs for SIR_dx_dt ODE
def make_inputs_helper(
        S0=50000-np.array([1]*2),
        I0=np.array([1]*2),
        R0=np.array([0]*2),
        CI0=np.array([0]*2),
        V0=np.array([0]*2),
        R0s=np.array([2,1.5]),
        alpha=0,
        gamma=1.0/3.0,
        alpha_hs_mult=[0]
        ):

    if not all(len(S0) == item for
               item in [len(I0),len(R0),len(CI0),len(V0),len(R0s)]):
        print "Inputs for state variables are not of the same length"
        return

    nlocs=S0.size

    betas=gamma*R0s
    initial_state=np.concatenate((S0,I0,R0,V0,CI0))
    if len(alpha_hs_mult) == 1:
        contact_mat=np.zeros((nlocs,nlocs)) + alpha/(nlocs-1.0)
        np.fill_diagonal(contact_mat, 1 - alpha)
    else:
        contact_mat = alpha_hs_mult

    return [initial_state,contact_mat,betas,gamma]

def make_inputs_helper_VE(
        S0=50000-np.array([1]*2),
        SV0=np.array([0]*2),
        I0=np.array([1]*2),
        R0=np.array([0]*2),
        CI0=np.array([0]*2),
        V0=np.array([0]*2),
        R0s=np.array([2,1.5]),
        alpha=0,
        theta=0.866,
        gamma=1.0/3.0,
        alpha_hs_mult=[0]
        ):

    if not all(len(S0) == item for
               item in [len(I0),len(SV0),len(R0),len(CI0),len(V0),len(R0s)]):
        print "Inputs for state variables are not of the same length"
        return

    nlocs=S0.size

    betas=gamma*R0s
    initial_state=np.concatenate((S0,SV0,I0,R0,V0,CI0))
    if len(alpha_hs_mult) == 1:
        contact_mat=np.zeros((nlocs,nlocs)) + alpha/(nlocs-1.0)
        np.fill_diagonal(contact_mat, 1 - alpha)
    else:
        contact_mat = alpha_hs_mult

    return [initial_state,contact_mat,betas,theta,gamma]



def SIR_dx_dt(state,t,contact_mat,betas,gamma,vac_start,
              vac_end,vac_daily_rates):
    '''Defining SIR System of Equations'''

    #print "SIR Time: %d" % (t)

    nlocs = state.size/5 # 4 states plus cumulative incidence
    S = np.array(state[0:nlocs])
    I = np.array(state[(nlocs):(2*nlocs)])
    R = np.array(state[(2*nlocs):(3*nlocs)])
    V = np.array(state[(3*nlocs):(4*nlocs)])
    N = S+I+R
    #print "S=%s,I=%s,R=%s,V=%s" % (S,I,R,V)

    ## make sure these don't go negative
    #I[I<0]=0.0
    #S[S<0]=0.0

    ## calculate the number of new infections
    new_inf = np.zeros(nlocs)
    for x in range(0,contact_mat.shape[0]):
        #        pdb.set_trace()
        new_inf[x] = np.sum((S[x]* contact_mat[x,:] * betas) * contact_mat[:,x].dot(I) / contact_mat[:,x].dot(N))

        #print "time is %s vac goes from %s to %s " % (t,vac_start,vac_end)
    if t >= vac_start and t < vac_end:
        #print "vaccinating at %s a total of %s" % (t,vac_daily_rates)
        new_vac = vac_daily_rates
    else:
        new_vac = np.zeros(vac_daily_rates.size)

    dS = -new_inf - (S/(N-V))*new_vac
    dI = new_inf - gamma*I
    dR = gamma*I + (S/(N-V))*new_vac
    dV = new_vac
    dCI = new_inf

    # print "dS=%s,dI=%s,dR=%s,dV=%s,CI=%s" % (dS,dI,dR,dV,dCI)
    #    pdb.set_trace()

    return np.concatenate((dS,dI,dR,dV,dCI))


def SIR_VE_dx_dt(state,t,contact_mat,betas,gamma,theta,vac_start,
              vac_end,vac_daily_rates):
    '''Defining SIR System of Equations'''
    ''' this one is for an imperfect all or nothing vaccine'''

    #print "SIR Time: %d" % (t)

    nlocs = state.size/6 # 5 states plus cumulative incidence
    S = np.array(state[0:nlocs])
    SV = np.array(state[(nlocs):(2*nlocs)])
    I = np.array(state[(2*nlocs):(3*nlocs)])
    R = np.array(state[(3*nlocs):(4*nlocs)])
    V = np.array(state[(4*nlocs):(5*nlocs)])
    N = S + I + R + SV
    #print "S=%s,I=%s,R=%s,V=%s" % (S,I,R,V)

    ## make sure these don't go negative
    #I[I<0]=0.0
    #S[S<0]=0.0

    ## calculate the number of new infections
    ## both in non-vaccinated suscpetiles and
    ## vaccinated susceptibles
    new_inf = np.zeros(nlocs)
    new_infV = np.zeros(nlocs)

    for x in range(0,contact_mat.shape[0]):
        new_inf[x] = np.sum((S[x]* contact_mat[x,:] * betas) * contact_mat[:,x].dot(I) / contact_mat[:,x].dot(N))
        new_infV[x] = np.sum((SV[x]* contact_mat[x,:] * betas) * contact_mat[:,x].dot(I) / contact_mat[:,x].dot(N))

        #print "time is %s vac goes from %s to %s " % (t,vac_start,vac_end)
    if t >= vac_start and t < vac_end:
        #print "vaccinating at %s a total of %s" % (t,vac_daily_rates)
        new_vac = vac_daily_rates
    else:
        new_vac = np.zeros(vac_daily_rates.size)

    dS = -new_inf - (S/(N-V))*new_vac
    dSV= -new_infV + (1-theta)*(S/(N-V))*new_vac
    dI =  new_inf + new_infV - gamma*I
    dR = gamma*I + theta*(S/(N-V))*new_vac
    dV = new_vac
    dCI = new_inf + new_infV

    # print "dS=%s,dI=%s,dR=%s,dV=%s,CI=%s" % (dS,dI,dR,dV,dCI)
    #    pdb.set_trace()

    return np.concatenate((dS,dSV,dI,dR,dV,dCI))


def get_final_attack_size(tc,n_states=5):
    ''' Gets the final attack size for all locations given output from SIR.dx.dt'''
    cols = tc.shape[1]
    nlocs = cols/n_states
    fs = [tc[-1:,x] for x in range(cols-nlocs,cols,1)]
    return fs

def get_ci_columns(tc,n_states=5):
    ''' Gets the final attack size for all locations given output from SIR.dx.dt'''
    cols = tc.shape[1]
    nlocs = cols/n_states
    cis = tc[:,range(cols-nlocs,cols,1)]
    return cis

# similar to expand.grid() in R
def expandgrid(*itrs):
   product = list(itertools.product(*itrs))
   return {'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))}

def two_pop_finalsize(individual,**options):
    """evalutates final size for a vaccine scenario:
    input:
    1. individual = vaccines in hotspot (will impute the other loc)
    2. options = all the deets for the simulation
    returns:
    tuple of the total funal size for the controlled epidemic
    [NOTE: this is used for optimization purposes]
    """

    my_alpha = options.get("alpha")
    my_N = options.get("N")
    my_nlocs = options.get("nlocs")
    R0 = options.get("R0")
    R1 = options.get("R1")
    # R2 = options.get("R2")
    my_vac_start = options.get("vac_start")
    my_times = options.get("times")
    my_hmax = options.get("hmax")
    vac_max = options.get("vac_max")

    individual = np.array([individual,vac_max-individual])

    # pdb.set_trace()

    ins = vf.make_inputs(my_N,
                         my_nlocs,
                         np.array([R0,R1]),
                         my_alpha,1.0/3,[0.0])

    # make vaccine allocation vector
    vac_daily = np.array(individual)

    tmp_run = spi.odeint(vf.SIR_dx_dt,
                         ins[0],
                         t=my_times,
                         hmax=my_hmax,
                         args=(ins[1],
                               ins[2],
                               ins[3],
                               my_vac_start,
                               my_vac_start+1,
                               vac_daily))

    final_size = sum(vf.get_final_attack_size(tmp_run))[0]

    return final_size

def two_pop_finalsize_helper(individual,alpha,N,nlocs,
             R0,R1,
             vac_start,
             times,
             hmax,
             vac_max):

    """helper evaluate function with all the defaults loaded """

    rc = two_pop_finalsize(individual,
             alpha=alpha,
             N=N,
             nlocs=nlocs,
             R0=R0,
             R1=R1,
             vac_start=vac_start,
             times=times,
             hmax=hmax,
             vac_max=vac_max)

    return rc
