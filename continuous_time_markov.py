#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Title: 
## Description: 
## Author:Guillaume Gay<elagachado AT  gmail DOT com>
## Commentary: Python version of X. Bressaud matlab code


 
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



from numpy import zeros, sign, floor, array, r_
from scipy.linalg import eig 

from pylab import figure, plot



### Initialisation des paramètres
k_a = 0.06
beta = 1
N = 10

def calc_rates(k_a, beta, N):
 
    Q=zeros((N + 1,) * 8)

    # On remplit les cases. Pour chaque etat, on calcule a quel taux on peut
    # passer aux differents etat qu'on peut atteindre en une etape. 

    for NAG in range(N + 1):
        for NAD in range(N + 1):
            for NBG in range(N + 1 ):
                for NBD in range(N + 1):
                    if not(NAD + NAG > N or NBG + NBD > N) :

                        # Calcul du taux de detachement.
                        forceA = NAD - NAG
                        forceB = NBD - NBG
                        if sign( forceA * forceB ) <= 0: 
                            M = abs(forceA - forceB)
                        else : 
                            M = - 4 / 3. # pour que 1+3M/8 = 1/2

                        k_d = k_a / (3 * ( 1 + 3 * M / 8. ))

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

    ### Transformation du tableau en une matrice 
    # Correspondance des indices. 

    compteur = 0
    AG = []
    AD = []
    BG = []
    BD = []
    for NAG in range(N + 1) :
        for NAD in range(N + 1) :
            for NBG in range(N + 1) :
                for NBD in range(N + 1) :
                    if not(NAD + NAG > N or NBG + NBD > N) :
                        AG.append(NAG)
                        AD.append(NAD)
                        BG.append(NBG)
                        BD.append(NBD)                    
                        compteur += 1

    AG = array(AG)
    AD = array(AD)
    BG = array(BG)
    BD = array(BD)

    taille=compteur   # taille de la matrice obtenue. 
    QQ=zeros((taille,taille)) # Initialisation 


    # La matrice 
    for i in range(taille):
        for j in range(taille):
            QQ[i,j] = Q[AG[i],AD[i],BG[i],BD[i],
                        AG[j],AD[j],BG[j],BD[j]]

    # Coefficients diagonaux

    for i in range(taille):
        QQ[i,i] = - sum(QQ[i,:])   
 

### Affichage des taux de transition partant de chaque état
# Bien sûr completement facultatif

 
# for i=1:compteur
#     ou=find(QQ(i,:));
#     alors=[AG(i)-1 AD(i)-1 BG(i)-1 BD(i)-1]
#     donc=[];
#     for j=1:length(ou)
#     donc=[donc;AG(ou(j))-1 AD(ou(j))-1 BG(ou(j))-1 BD(ou(j))-1 QQ(i,ou(j))];
#     end
#     donc
#     #Alors= [AG(i)+AD(i) BG(i)+BD(i)]
#     #Donc= [donc(:,1)+donc(:,2) donc(:,3)+donc(:,4)]
# end

### Recherche du noyau, i.e. de la mesure invariante
 
#[V,D] = eig(QQ'); # V les vecteurs propres, D la matrice diagonale. 
#en python:

valeurs_propres, vecteurs_propres = eig(QQ.T)
 
# Attention, QQ.T est la transposee de QQ. 
# En effet, on a toujours QQ.1 = 0 (vu ce qu'on a mis sur la diagonale)
# et nous, on cherche le vecteur V tel que V'. QQ=0, soit QQ'.V=0
 
# Determination de l'indice de la valeur propre qui est nulle

 
for i, v in enumerate(valeurs_propres):
    if abs(v)< 1e-10:
         uu = i
 
# dot(QQ, ones(taille)) # pour verifier que la somme des lignes est nulle
# dot(QQ.T, V[:,uu]) # pour verifier que le vecteur obtenu est bien le noyau
 
# Le vecteur ppre normalise donne la proba invariante

mes = vecteurs_propres[:,uu]/sum(vecteurs_propres[:,uu]) 
 
# Pour voir mieux, je regroupe les etats en termes de la distance d
# Pour l'instant, je classe violemment en fonction de la difference des
# forces. Mais il faudra affiner. 

 
M = floor(abs( AG + BD - AD - BG ) )
 
# Visualisation des configurations associée a chaque valeur de M, et de la
# proba invariante de cet configuration
 
#Voir=[mes, M.T, AG.T, AD.T, BG.T, BD.T]

 
# On cumule les proba des etats donnant la meme valeur à M
pp = zeros(2*N + 1)

for i, m in enumerate(mes):
    pp[M[i]] += m
 
# Affichage du resultat 

#Distance=[pp (0:2*N)']

 
# Tracé de la loi obtenue: 
# figure(1)
# plot(r_[0:2*N+1], pp, '*-')
