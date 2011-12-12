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
beta = 0.5
N = 4
kappa = 3/8.
d_alpha = 1/3.

betas = linspace(0,1,10) 
d_alphas = logspace(-2, 2, 10)


def explore(betas = betas, d_alphas = d_alphas):
    
    espace = StateSpace(N)
    measures = []
    for beta in betas:
        for d_alpha in d_alphas:
            print 'beta = %.3f' %beta
            generateur = ChainGenerator(espace, beta = beta,
                                        d_alpha = d_alpha)
            measures.append(real(generateur.measure))
    return measures

def all_defects_probas(espace, measures,
                       d_alphas = d_alphas,
                       betas = betas, N = N):

    all_defect_probs = {'amphiteliques':zeros(len(measures)),
                        'monoteliques':zeros(len(measures)),
                        'meroteliques':zeros(len(measures)),
                        'synteliques':zeros(len(measures)),
                        'unattached':zeros(len(measures))}
    all_p_cor, all_p_err, all_p_force = (zeros((len(measures),
                                                2 * N + 1)),) * 3
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
            xscale('log')
            axis((d_alphas.min(), d_alphas.max(),
                  betas.min(), betas.max()))
            if i == 1 or i == 3:
                ylabel(r'Orientation effect $\beta$')
            if i == 3 or i == 4:
                xlabel(r'$d_\alpha$')

            if i == 2:
                xticks(())
                yticks(())
            if i == 4:
                yticks(())
            if i == 1:
                xticks(())

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



class ChainGenerator():

    def __init__(self, state_space, k_a = k_a, beta = beta,
                 N = N, kappa = kappa, d_alpha = d_alpha):
        '''
        Intanciation method for a ChainGenerator
        The following attributes are caclulated :

        self.Q : a 8 * taille array (vector in / vector out)
        self.QQ: the actual matrix which is the generator of the chain
        self.vecteurs_prore, self.valeurs_propres: result of eig(self.QQ.T)
        self.measure: the invariant measure
        self.show_invariant(display = True) 
        '''

        valid_states = state_space.valid_states
        self.calc_rates(valid_states, k_a, beta,
                        N, kappa, d_alpha)
        self.get_matrix(valid_states)
        self.find_kernel()

    
    def calc_rates(self, valid_states, k_a,
                   beta, N, kappa, d_alpha):
 
        self.Q=zeros((N + 1,) * 8)

        # On remplit les cases. Pour chaque etat, on calcule a quel taux on peut
        # passer aux differents etat qu'on peut atteindre en une etape. 

        for psi in valid_states.T : #On itère sur tous les etats valides

            NAG, NAD, NBG, NBD = psi
            if NAG + NBD >=  NAD +  NBG:
                sigma = 1
            else:
                sigma = -1

            d_eq = 1 + sigma * 3. * (NBD - NBG + NAG - NAD) / (2. * N)
            # Calcul du taux de detachement.
            k_d = k_a * d_alpha / d_eq
            
            #la vieille definition
            # forceA=NAD-NAG
            # forceB=NBD-NBG
            # if sign(forceA*forceB)<=0 :
            #     M=abs(NAG+NBD - NAD - NBG); #i.e. somme des valabs des forces
            # else:
            #     M=-1/3 # pour que 1+3M/8 = 1/2
            # k_d= k_a / (3. * ( 1 + 3. * M / 8))

            # Taux des differents detachements
            if NAG > 0:
                self.Q[NAG, NAD, NBG, NBD,
                       NAG - 1, NAD, NBG, NBD] = k_d * NAG

            if NAD > 0:
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD - 1, NBG, NBD] = k_d * NAD

            if NBG > 0:
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD, NBG - 1, NBD] = k_d * NBG

            if NBD > 0:
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD, NBG, NBD - 1] = k_d * NBD

            # Taux d'attachement (avec proba pour droit ou gauche)
            if NAG + NAD == 0 :  # sinon division par zero
                Pe = 0.5
            else :
                Pe = 0.5 + beta * (NAG - NAD) / (2. * (NAG + NAD)) 

            if NAG < N :
                self.Q[NAG, NAD, NBG, NBD,
                       NAG + 1, NAD, NBG, NBD] = k_a * (N - NAG - NAD) * Pe
            if NAD < N :
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD + 1, NBG, NBD] = k_a * (N - NAG - NAD
                                                        ) * ( 1 - Pe )

            if NBG + NBD == 0 : # sinon division par zero
                Pe = 0.5
            else : 
                Pe = 0.5 + beta * (NBG - NBD) / (2. * (NBG + NBD))
            if NBG < N :
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD, NBG + 1, NBD] = k_a * (N - NBG - NBD
                                                        ) * Pe 
            if NBD < N :
                self.Q[NAG, NAD, NBG, NBD,
                       NAG, NAD, NBG, NBD + 1] = k_a * (N - NBG - NBD
                                                        ) * (1 - Pe )

    def get_matrix(self, valid_states):

        ### Transformation du tableau en une matrice 
        taille = valid_states.shape[1]
        self.QQ=zeros((taille,taille)) # Initialisation 
        # La matrice 
        for i,j in ndindex((taille, taille)):
            self.QQ[i,j] =self.Q[valid_states[0,i],
                                 valid_states[1,i],
                                 valid_states[2,i],
                                 valid_states[3,i],
                                 valid_states[0,j],
                                 valid_states[1,j],
                                 valid_states[2,j],
                                 valid_states[3,j]]
        # Coefficients diagonaux
        for i in range(taille):
            self.QQ[i,i] = - sum(self.QQ[i,:])   

    def find_kernel(self):

        ### Recherche du noyau, i.e. de la mesure invariante
        self.valeurs_propres, self.vecteurs_propres = eig(self.QQ.T)

        # Attention, QQ.T est la transposee de QQ. 
        # En effet, on a toujours QQ.1 = 0 (vu ce qu'on a mis sur la diagonale)
        # et nous, on cherche le vecteur V tel que V'. QQ=0, soit QQ'.V=0

        # Determination de l'indice de la valeur propre qui est nulle
        for i, v in enumerate(self.valeurs_propres):
            if abs(v)< 1e-10:
                 uu = i
        # dot(QQ, ones(taille)) # pour verifier que la somme des lignes
        # est nulle
        # dot(QQ.T, V[:,uu]) # pour verifier que le vecteur obtenu
        # est bien le noyau
        # Le vecteur ppre normalise donne la proba invariante
        self.measure = self.vecteurs_propres[:,uu]/sum(
            self.vecteurs_propres[:,uu]) 
    
def show_invariant(measure, espace, display = True):

    p_cor = zeros(2*N + 1)
    p_err = zeros(2*N + 1)
    p_force = zeros(2*N + 1)

    defect_probs = {}

    for key, defect in espace.defects.items():
        defect_probs[key] = sum(measure * defect)


    # On cumule les proba des etats donnant les memes valeurs
    for i, m in enumerate(real(measure)):
        p_cor[espace.sum_corrects[i]] += m
        p_err[espace.sum_errones[i]] += m
        p_force[espace.force_totale[i]] += m

    if display == True:
        for key, proba in defect_probs.items():
            print 'P(%s) = %.03f' %(key, proba)
        #Trace de la loi obtenue: 
        figure(1)
        plot(r_[0:2*N+1], p_cor, 'gv-', label='Corectly attached')
        plot(r_[0:2*N+1], p_err, 'ro-', label='Erroneously attached')
        figure(2)
        plot(r_[0:2*N+1], p_force, 'ko-')
        title('Force totale')
    probas = defect_probs, p_cor, p_err, p_force
    return probas

    
class StateSpace():
    '''
    For a given number of attachment site per chromosomes, calculates the
    valid state space (such that NAD + NAG <= N et NBG + NBD <= N)

    attributes:
    self.valid_states a (N, number_of_states) array ({\Psi} in the text)
    self.new_states: a (N, number_of_states) array ({\Phi} in the text)
    self.defects = {'amphitelics':amphiteliques,
                    'monotelics':monoteliques,
                    'merotelics':meroteliques,
                    'syntelics':synteliques,
                    'unattached':unattached}
    where amphiteliques, monoteliques, etc. are 1D array such that
    if amphitelique[i]==1, then the state self.valid_states[i] is amphitelic.
    '''

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


if __name__ == "__main__":

    print 'usage (for example):'
    print 'espace = StateSpace(N) #N is the number of attachment sites'
    print 'generateur = ChainGenerator(espace, beta = beta,'
    print '                            d_alpha = d_alpha)'
    print 'generateur.show_invariant(espace, display = True)'
    print 'see also the explore() function'
    
