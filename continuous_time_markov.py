#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:Guillaume Gay<elagachado AT  gmail DOT com>
## Commentary: Python version of X. Bressaud matlab code


from numpy import zeros, ones, sign, floor, array, r_, linspace
from numpy import ndindex, hstack, identity, real, logspace
from scipy.linalg import eig 

from pylab import *
 
#####################################
### Modele
 
# chaque site d'attachement peut etre dans 3 etat rien, gauche, droite
# Il y a Mk sites d'attachement sur le centromere A et Mk sur le B.
 
# L etat du systeme est decrit par
# Le nombre de site d'attachements de A dans letat gauche NAG
# Le nombre de site d'attachements de A dans letat droite NAD
# Le nombre de site d'attachements de B dans letat gauche NBG
# Le nombre de site d'attachements de B dans letat droite NBD
  
 
# On veut evaluer les taux de transition en fonction de l'etat.
# Pour chacune des composantes, le raisonnement est le même:

 
# Regardons NAG

 
# NAG peut diminuer de 1. 
 
# Cela arrive avec un taux proportionnel a NAG
# C'est aussi fonction de la "distance". 
# On ecrit: k_d = k_a d_alpha/ d
# ou d_alpha de l'ordre de 3 d_0 (d_0 a lequilibre)
# et d la distance réelle entre les deux centromeres.
# On ecrit d=(1 + 3M/8) d_0, 
# de maniere a ce que si M=0 on obtienne d_0 et si M=8, 4 d_0
# Pour evaluer M on distingue les cas ou les forces exercées sur A et B
# sont de même signe, de signes opposés 
# Si elles sont de signe opposé, on pose M = somme des val abs de ces forces
# cad: M= abs(NAG + NBD - NAD - NBG)
# Sinon on pose M=-4/3 pour avoir d= d_0/2  
# Si les deux sont nulles, on s'en fout puisque personne attaché. 
 
# NAD+NAG peut augmenter de 1.
 
# Le taux est proportionnel au nombre de site d'attachements libres,
# i.e. a N-NAG-NAD
# La constante est k_a
# Mais alors c'est NAG ou NAD qui augmente de 1 ; la proba que ce soit NAG
# est donnée par 1/2 + beta * (NAG - NBG) / 2 (NAG + NAD)
#####################################

### Calcul des taux
# On definit un tableau Q ecrit en dimension 8 ; mais il faut voir comme
# une matrice dont les lignes et les colonnes sont indexees par des
# quadruplets decrivant les etats possibles du systeme. 

 
# Q(NAG,NAD,NBG,NBD,NAG',NAD',NBG',NBD') est le taux auquel le systeme
# passe de l'etat (NAG,NAD,NBG,NBD) a l'etat (NAG',NAD',NBG',NBD').  
 
# Remarque: or(NAD-1+NAG-1>N,NBG-1+NBD-1>N)==0 permet de ne considere que
# des etats raisonnables, i.e. ou le nombre de site d'attachements utilisés ne
# depasse pas N





### Initialisation des paramètres
k_a = 0.06
beta = 1.
N = 4
kappa = 3/8.
d_alpha = 1/3.

betas = linspace(0,1,20) 
d_alphas = logspace(-2, 1, 20)


def explore(betas = betas, d_alphas = d_alphas):

    espace = StateSpace(N)
    pcs, pes, measures = [], [], []

    for d_alpha in d_alphas:
        print 'd_alpha = %.3f' %d_alpha
        for beta in betas:
            mes = find_invariant(espace.valid_states,
                                 d_alpha = d_alpha, beta = beta)
            measures.append(mes)

    return measures




def all_defects(espace, measures, d_alphas = d_alphas, betas = betas, N = N):

    all_defect_probs = {'amphiteliques':zeros(len(measures)),
                        'monoteliques':zeros(len(measures)),
                        'meroteliques':zeros(len(measures)),
                        'synteliques':zeros(len(measures)),
                        'unattached':zeros(len(measures))}

    all_p_cor, all_p_err, all_p_force = (zeros((len(measures), 2 * N + 1)),) * 3

    for i, mes in enumerate(measures):
        probas = show_invariant(mes, espace, display = False)
        defect_probs, p_cor, p_err, p_force = probas

        for key, prob in defect_probs.items():
            all_defect_probs[key][i] = prob
        all_p_cor[i, :] = p_cor 
        all_p_err[i, :] = p_err 
        all_p_force[i, :] = p_force 

    for key, probs in all_defect_probs.items():
        all_defect_probs[key] = probs.reshape((d_alphas.size, betas.size))

    all_p_cor = all_p_cor.reshape((d_alphas.size, betas.size, 2 * N + 1))
    all_p_err = all_p_err.reshape((d_alphas.size, betas.size, 2 * N + 1))
    all_p_force = all_p_force.reshape((d_alphas.size, betas.size, 2 * N + 1))
    
    probas = all_defect_probs, all_p_cor, all_p_err, all_p_force
    return probas

def show_all_defects(probas, d_alphas = d_alphas, betas = betas):

    defect_probs, p_cor, p_err, p_force = probas
    
    total = zeros((d_alphas.size, betas.size))

    nd = zeros(d_alphas.size + 1)
    nb = zeros(betas.size + 1)

    nd[:d_alphas.size] = d_alphas
    nd[-1] = 2 * d_alphas[-1] - d_alphas[-2]
    nb[:betas.size] = betas
    nb[-1] = 2 * betas[-1] - betas[-2]
    d_alphas, betas = nd, nb

    figure()
    i = 1
    for key, defect in defect_probs.items():
        if 'amphi' not in key:
            subplot(2,2,i)
            total += defect
            pcolor(d_alphas, betas, defect,
                   vmin = 0., vmax = 1.)
            colorbar()
            xscale('log')
            ylabel(r'Orientation effect $\beta$')
            xlabel(r'$d_\alpha$')
            axis((d_alphas.min(), d_alphas.max(), 0, 1))
            title(key)
            i+=1

    figure()
    pcolor(d_alphas, betas, total)
    colorbar()
    xscale('log')
    axis((d_alphas.min(), d_alphas.max(), 0, 1))
    ylabel(r'Orientation effect $\beta$')
    xlabel(r'$d_\alpha$')
    title('Total')



def find_invariant(valid_states, k_a = k_a, beta = beta,
                   N = N, kappa = kappa, d_alpha = d_alpha):
    
    Q = calc_rates(valid_states, k_a, beta, N, kappa, d_alpha)
    QQ = get_matrix(valid_states, Q)
    mes = find_kernel(QQ)
    return real(mes)

def calc_rates(valid_states, k_a, beta, N, kappa, d_alpha):
 
    Q=zeros((N + 1,) * 8)

    # On remplit les cases. Pour chaque etat, on calcule a quel taux on peut
    # passer aux differents etat qu'on peut atteindre en une etape. 

    for psi in valid_states.T : #On itère sur tous les etats valides

        NAG, NAD, NBG, NBD = psi
        # Calcul du taux de detachement.
        forceA = NAD - NAG
        forceB = NBD - NBG
        #With Fk as unit force and d_0 the unit distance:
        if sign( forceA * forceB ) <= 0:
            d_eq =  1. + abs(forceA - forceB) / kappa 
            k_d = k_a * d_alpha / d_eq
        else:
            k_d = k_a * 5. 

        # Taux des differents detachements
        if NAG > 0:
            Q[NAG, NAD, NBG, NBD,
              NAG - 1, NAD, NBG, NBD] = k_d * NAG

        if NAD > 0:
            Q[NAG, NAD, NBG, NBD,
              NAG, NAD - 1, NBG, NBD] = k_d * NAD

        if NBG > 0:
            Q[NAG, NAD, NBG, NBD,
              NAG, NAD, NBG - 1, NBD] = k_d * NBG

        if NBD > 0:
            Q[NAG, NAD, NBG, NBD,
              NAG, NAD, NBG, NBD - 1] = k_d * NBD
            
        # Taux d'attachement (avec proba pour droit ou gauche)
        if NAG + NAD == 0 :  # sinon division par zero
            Pe = 0.5
        else :
            Pe = 0.5 + beta * (NAG - NAD) / (2. * (NAG + NAD)) 

        if NAG < N :
            Q[NAG, NAD, NBG, NBD,
              NAG + 1, NAD, NBG, NBD] = k_a * (N - NAG - NAD) * Pe
        if NAD < N :
            Q[NAG, NAD, NBG, NBD,
              NAG, NAD + 1, NBG, NBD] = k_a * (N - NAG - NAD) * (1-Pe)

        if NBG + NBD == 0 : # sinon division par zero
            Pe = 0.5
        else : 
            Pe = 0.5 + beta * (NBG - NBD) / (2. * (NBG + NBD))
        if NBG < N :
            Q[NAG, NAD, NBG, NBD,
              NAG, NAD, NBG + 1, NBD] = k_a * (N - NBG - NBD) * Pe 
        if NBD < N :
            Q[NAG, NAD, NBG, NBD,
                  NAG, NAD, NBG, NBD + 1] = k_a * (N - NBG - NBD) * (1 - Pe )

    return Q


def get_matrix(valid_states, Q):

    ### Transformation du tableau en une matrice 
    taille = valid_states.shape[1]
    QQ=zeros((taille,taille)) # Initialisation 
    # La matrice 
    for i in range(taille):
        for j in range(taille):
            indices = tuple(hstack((valid_states[:,i],
                                    valid_states[:,j])))
            QQ[i,j] = Q[indices]
            
    # Coefficients diagonaux
    for i in range(taille):
        QQ[i,i] = - sum(QQ[i,:])   

    return QQ

def find_kernel(QQ):

    ### Recherche du noyau, i.e. de la mesure invariante
    valeurs_propres, vecteurs_propres = eig(QQ.T)

    # Attention, QQ.T est la transposee de QQ. 
    # En effet, on a toujours QQ.1 = 0 (vu ce qu'on a mis sur la diagonale)
    # et nous, on cherche le vecteur V tel que V'. QQ=0, soit QQ'.V=0

    # Determination de l'indice de la valeur propre qui est nulle
    for i, v in enumerate(valeurs_propres):
        if abs(v)< 1e-10:
             uu = i
    # dot(QQ, ones(taille)) # pour verifier que la somme des lignes est nulle
    # dot(QQ.T, V[:,uu]) # pour verifier que le vecteur obtenu est bien le noyau        # Le vecteur ppre normalise donne la proba invariante
    mes = vecteurs_propres[:,uu]/sum(vecteurs_propres[:,uu]) 
    return mes
    
def show_invariant(mes, espace, display = True):
    
    p_cor = zeros(2*N + 1)
    p_err = zeros(2*N + 1)
    p_force = zeros(2*N + 1)
    defect_probs = {}
    
    for key, defect in espace.defects.items():
        defect_probs[key] = sum(mes * defect)


    # On cumule les proba des etats donnant les memes valeurs
    for i, m in enumerate(real(mes)):
        p_cor[espace.sum_corrects[i]] += m
        p_err[espace.sum_errones[i]] += m
        p_force[espace.force_totale[i]] += m

    if display == True:
        for key, proba in defect_probs.items():
            print 'P(%s) = %.03f' %(key, proba)

        #Trace de la loi obtenue: 
        figure(1)
        plot(r_[0:2*N+1], p_cor, 'gv-')
        plot(r_[0:2*N+1], p_err, 'ro-')
        figure(2)
        plot(r_[0:2*N+1], p_force, 'ko-')


    probas = defect_probs, p_cor, p_err, p_force
    return probas

    
class StateSpace():


    def __init__(self, N):

        self.N = N
        self.get_valid_states()
        self.get_correct_states()
        self.get_defects()
        
    def get_valid_states(self):

        N = self.N

        #On commence par trouver les combinaisons d'etats valides
        #donc telles que NAD + NAG <= N et NBG + NBD <= N)
        compteur = 0
        AG = []
        AD = []
        BG = []
        BD = []

        for NAG, NAD, NBG, NBD in ndindex(N+1, N+1, N+1, N+1):
            if not(NAD + NAG > N or NBG + NBD > N) :
                AG.append(NAG)
                AD.append(NAD)
                BG.append(NBG)
                BD.append(NBD)                    

        AG = array(AG)
        AD = array(AD)
        BG = array(BG)
        BD = array(BD)
        self.valid_states = array([AG, AD, BG, BD])


    def get_correct_states(self):

        AG, AD, BG, BD = self.valid_states
        AC, AE, BE, BC = self.valid_states.copy()

        switch = AG + BD < AD + BG
        
        AC[switch] = AD[switch]
        AE[switch] = AG[switch]

        BC[switch] = BG[switch]
        BE[switch] = BD[switch]
        self.new_states = array([AC, AE, BE, BC])


    def get_defects(self):

        AC, AE, BE, BC = self.new_states

        unattached = zeros(AE.shape)
        unattached[AC +AE + BE + BC == 0] += 1

        meroteliques = zeros(AE.shape)
        meroteliques[AC * AE != 0 ] += 1
        meroteliques[BC * BE != 0 ] += 1
        # On compte les doubles comme un seul..
        meroteliques[meroteliques == 2] = 1


        monoteliques = zeros(AE.shape)
        monoteliques[AC + AE == 0] += 1
        monoteliques[BC + BE == 0] += 1
        monoteliques[AC + AE + BC + BE == 0] = 0
        monoteliques[meroteliques != 0] = 0

        synteliques = zeros(AE.shape)
        synteliques[AC * BE != 0] += 1
        synteliques[AE * BC != 0] += 1
        synteliques[meroteliques != 0] = 0
        
        amphiteliques = ones(AE.shape)
        amphiteliques[AE + BE != 0] = 0
        amphiteliques[AC * BC == 0] = 0


        self.defects = {'amphiteliques':amphiteliques,
                        'monoteliques':monoteliques,
                        'meroteliques':meroteliques,
                        'synteliques':synteliques,
                        'unattached':unattached}

        self.sum_corrects = AC + BC
        self.sum_errones = AE + BE
        self.force_totale = abs(AC + BC - AE - BE )
