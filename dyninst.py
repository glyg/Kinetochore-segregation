#!/usr/bin/python
# -*- coding: utf-8 -*-


### Adapted from Civelekoglu-Scholey et al. Biophys.J 90(11) 2006
### doi: 10.1529/biophysj.105.078691

"""
This module calculates the transition matrix of attachment/detachment mechanism
In its present form, it is limited to Mk = 4 and merotelic attachment is not
accounted for.



"""




from numpy import array, logspace, dot, mean, linspace, log, ndindex, zeros, exp
from scipy import linalg, column_stack, comb
from scipy.interpolate import splrep, splev


__all__ = ["mean_attachment", "delta_Mk", "explore_delta", "stat_dist"]


def _prob(freq, dt):
    return 1 - exp( - freq * dt )

def _prob_4j(Pd, j):

    return float(comb(4,j, exact=1)) * (Pd**(4 - j)) * (1 - Pd)**j

def _prob_0j(Pa,j):

    return _prob_4j(1-Pa, j)

def _prob_j4(Pa, Pd, j):

    return Pa**(4 - j) * (1 - Pd)**j

def _prob_j0(Pa, Pd, j):

    return _prob_j4(1-Pa, 1-Pd, j)

def _prob_33(Pa, Pd):

    return (1 - Pa) * (1 - Pd)**3 + 3 * Pa * Pd * (1 - Pd)**2

def _prob_11(Pa, Pd):

    return _prob_33(Pd, Pa)

def _prob_32(Pa, Pd):

    return 3 * (1 - Pa) * Pd * (1 - Pd)**2 + 3 * Pa * Pd**2 * (1 - Pd)

def _prob_12(Pa, Pd):

    return _prob_32(Pd, Pa)

def _prob_31(Pa, Pd):

    return 3 * (1 - Pa) * Pd**2 * (1 - Pd) + Pa * Pd**3 

def _prob_13(Pa, Pd):

    return _prob_31(Pd, Pa)

def _prob_23(Pa, Pd):

    return 2 * (1 - Pa) * Pa * (1 - Pd)**2 + 2 * Pa**2 * Pd * (1 - Pd)

def _prob_21(Pa, Pd):

    return _prob_23(Pd, Pa)

def _prob_22(Pa,Pd):

    return (1 - Pd)**2 * (1 - Pa)**2 + 4 * Pd * Pa * (1 - Pd) * (1 - Pa) + Pd**2 * Pa**2
    
def trans_matrix(Pa, Pd):
    """returns the transition matrix for probabilities Pa and Pd """
    T = zeros((5, 5), float)
    for j in range(5):
        T[4, j] = _prob_4j(Pd, j)
        T[j, 4] = _prob_j4(Pa, Pd, j)
        T[0, j] = _prob_0j(Pa, j)
        T[j, 0] = _prob_j0(Pa, Pd, j)

    T[1, 1] = _prob_11(Pa, Pd)
    T[2, 1] = _prob_21(Pa, Pd)
    T[3, 1] = _prob_31(Pa, Pd)
    T[1, 2] = _prob_12(Pa, Pd)
    T[2, 2] = _prob_22(Pa, Pd)
    T[3, 2] = _prob_32(Pa, Pd)
    T[1, 3] = _prob_13(Pa, Pd)
    T[2, 3] = _prob_23(Pa, Pd)
    T[3, 3] = _prob_33(Pa, Pd)

    return T


def stat_dist(Pa, Pd):
    """Returns the unitary eigen vector of the transition matrix
    """
    T = trans_matrix(Pa,Pd)
    (evals, evects) = linalg.eig(T.T)
    for i in range(5):
        if 1.0 - 1e-15 < evals[i] <= 1.0 +1e-15 :
            stat_dist =  abs(evects[:,i]) / sum(abs(evects[:,i]))
            break
        elif i == 4:
            print "unitary eigenvalue not found !"
            return 0

    check = max(abs(dot(stat_dist, T) - stat_dist))

    if check > 1e-15:
        print "please examine the vector"

    return stat_dist
    
    
def explorons():

    alphas = logspace(-2.,2.,base = 10.0)

    probs = []
    coefs = []
    means = []
    
    for alpha in alphas:

        Pa = _prob(0.01*alpha, 1.)
        Pd = _prob(0.01, 1.)
        probs.append(Pd / Pa)

        sd = stat_dist(Pa, Pd)
        if array(sd).any():
            coefs.append(sd)
            means.append(dot(sd,range(5)))
        else:
            print alpha
            coefs.append(zeros(5))
            means.append(0)
    probs = array(probs)
    coefs = array(coefs)
    means = array(means)
    i_mean = splrep(log(alphas), means)
    i_coefs = []
    for k in range(coefs.shape[1]):
        i_coefs.append(splrep(log(alphas), coefs[:,k]))
    return i_coefs,  i_mean

def mean_attachment(ratio):
    '''
    ratio is attachment over detachment frequencies

    returns the mean number of attached sites for the given ratio
    
    '''
    i_coefs, i_mean = explorons()

    mean_Mk = splev(log(ratio), i_mean)

    return mean_Mk

def delta_Mk(ratio):

    i_coefs, i_mean = explorons()
    
    delta = zeros(9)
    sigma = zeros(9)
    products = zeros((5,5))
    for i, j in ndindex((5,5)):
        products[i,j] = splev(log(ratio), i_coefs[i]) * splev(log(ratio), i_coefs[j])
        delta[i - j + 4] += products[i,j]
        sigma[i + j] += products[i,j]

    abs_delta = delta[-5:]*2.
    m_delta = sum(abs_delta*arange(5))

    m_sigma = sum(sigma*arange(9))

    return delta, sigma, m_delta, m_sigma
    
        
def explore_delta():

    alphas = logspace(-2.,2.,base = 10.0)

    sigmas = []
    deltas = []
    m_deltas = []
    m_sigmas = []
    for alpha in alphas:
        delta, sigma, m_delta, m_sigma = delta_Mk(alpha)
        sigmas.append(sigma)
        deltas.append(delta)
        m_sigmas.append(m_sigma)
        m_deltas.append(m_delta)

    
    return alphas, sigmas, deltas, m_deltas, m_sigmas

    
