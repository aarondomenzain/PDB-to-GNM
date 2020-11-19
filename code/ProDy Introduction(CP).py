#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 22:12:37 2020
@author: aarondomenzain
"""


from prody import *
import matplotlib.pyplot as plt
import numpy as np


#Compute free energy from given PDB structure
def gnmanalysis (protPDB, model_number):
    p = parsePDB(protPDB, model = model_number)
    
    # Alpha Carbons
    calphas = p.select('calpha')
    
    # P Atoms
    phospho = p.select(' name "P.*" ')
    
    # Checks if calphas and phospho instances are empty. If not, it builds up the whole network
    if phospho is not None and calphas is not None:
        network = calphas + phospho
    
    else:
        
        # Build network with P atoms only 
        if phospho is not None:
        
            # Add P atoms to network
            network = phospho

        # Build network with C alpha atoms only
        if calphas is not None:
            
            # Add Calpha atoms to network
            network = calphas
    
    #Gaussian Network Modelling
    gnm = GNM(protPDB)
    gnm.buildKirchhoff(network, cutoff = 10.0, gamma = 1.0)

    #Kirchoff Matrix
    K = gnm.getKirchhoff()

    #Calculate Normal Modes
    gnm.calcModes(1000) #All normal modes (default is 20)

    #Eigenvalues and eigenvectors
    Eigvals = gnm.getEigvals() #This is a vector containing eigenfrequencies
    Eigvecs = gnm.getEigvecs() #This is a matrix, each column is an eigenvector (normal mode)
    

    #Compute Free Energy in KT units (F)
    E = np.sum(np.log(Eigvals))
    
    return (E, Eigvals, Eigvecs)

#Iterates over the PDB structure number model 
def gnmanalysis_models(protname, n_models):
    #Initialize a list of zeros of length (n_models)
    listEnergies = np.zeros(n_models)
    
    for n in range (n_models):
        
        Energy, Eigvals, Eigvecs = gnmanalysis(protname, n+1)
        
        #No encontré manera bonita de iniciar matrices de ceros porque necesito conocer el número de eigenvectores y su tamaño
        if n ==0:
            Neigvals = len(Eigvals) 
            Neigvecs = len(Eigvecs)
            listEigvals = np.zeros([n_models, Neigvals])
            #Calculate dimensions of eigenvectors matrix
            [j, k] = Eigvecs.shape
            #Initialize zeros array
            listEigvecs = np.zeros([n_models, j, k])
        listEnergies[n] = Energy
        listEigvals[n] = Eigvals
        listEigvecs[n] = Eigvecs
        
    np.save('../results/' + protname + "_model_energies.npy", listEnergies)
    np.save('../results/' + protname + "_eigenvalues.npy", listEigvals)
    np.save('../results/' + protname + "_eigenvectors.npy", listEigvecs)
    # return (listEnergies, listEigvals, listEigvecs)

def hist_plot(frequencies, protname, color):
    
    # List of colors
    # 'lighblue'
    # 'lightgreen'
    # 'salmon'
    
    plt.hist(frequencies, density=1, color = color, label = protname)
    plt.xlabel("Eigenfrecuencias (U.A.)")
    plt.ylabel("Densidad de estados g(w)")
    plt.xlim(0, 25)
    plt.ylim(0, 0.1)
    plt.legend()
    plt.show()


#####################      INSTRUCTIONS       ##########################
    
# GNM Normal Mode Analysis
    

# PDB list to analyze:
    #Unbounded 
#     '1AAF' : Unbounded NC from HIV-1 (20 models) color = 'gray'
#     '1ESY' : Unbounded SL2 from HIV-1 (19 models) color = 'lightblue'
#     '1BNO' : Unbounded SL3 from HIV-1 (11 models) color = 'lightgreen'
    
    #Complexed
#     '1F6U' : NC+SL2 complex from HIV-1 (20 models) color = 'blue'
#     '1A1T' : NC+SL3 complex from HIV-1 (25 models) color = 'green'


# PDB_list = ['1AAF', '1ESY', '1BNO', '1F6U', '1A1T' ]
# N_list = [ 20, 19, 11, 20, 25 ]





# STEPS 
    
# 1. RUN THE SCRIPT AS IS TO LOAD FUNCTIONS


# 2. SET SYSTEM TO ANALYZE

    # protPDB = PDB to analyze (ex. prot = '1aaf')
    # N = number of structural models in PDB file  (ex. N = 20  -for 1aaf- )
    

# 3. EXECUTE 
    
  # Calculation of Gibbs free energy distribution of the PDB file, eigenfrequencies and eigenvectors
    # E, eigvals, eigvecs = gnmanalysis_models(prot, N)
    
    
 # Save data to Numpy files
    # np.save(protPDB + "_model_energies.npy", E)
    # np.save(protPDB + "_eigenvalues.npy", eigvals)
    # np.save(protPDB + "_eigenvectors.npy", eigvecs)

     # For loading Numpy files
        # E = np.load( protPDB + "_model_energies.npy" )
        # eigvals = np.load( protPDB + "_eigenvalues.npy" )
        # eigvecs = np.load( protPDB + "_eigvecs.npy" )

    
# #Quick histogram of first model
# plt.hist(eigvals[0], label= 'Distribución de estados', density= True, color= 'Green')
# plt.xlabel('Eigenfrecuencias (U.A.)')
# plt.ylabel('Densidad de estados')
# plt.ylim(0, 0.2)
# plt.grid(True)
# plt.legend()

# NCenergy = np.load("1aaf_model_energies.npy")
# SL2energy = np.load("1esy_model_energies.npy")
# SL3energy = np.load("1bn0_model_energies.npy")
# NCSL2energy = np.load("1f6u_model_energies.npy")
# NCSL3energy = np.load("1a1t_model_energies.npy")

# NCeigvals = np.load("1aaf_eigenvalues.npy")
# SL2eigvals = np.load("1esy_eigenvalues.npy")
# SL3eigvals = np.load("1bn0_eigenvalues.npy")
# NCSL2eigvals = np.load("1f6u_eigenvalues.npy")
# NCSL3eigvals = np.load("1a1t_eigenvalues.npy")




# # Load data 
# complex = [NCenergy, SL2energy, SL3energy, NCSL2energy, NCSL3energy]
# names = ['NC', 'SL2', 'SL3', 'NC+SL2', 'NC+SL3']

# highs = [NCenergy, NCSL2energy, NCSL3energy]ß
# highnames = ['NC', 'NC+SL2', 'NC+SL3']

# lows = [SL2energy, SL3energy]
# lownames = ["SL2", "SL3"]


# ## Free energy distributions
# plt.boxplot(complex, notch = False, showmeans = True, meanline = True)#, showfliers = False)
# plt.xticks([1, 2, 3, 4, 5], names)
# plt.ylabel('Energía libre (KT)')
# plt.show()
    
# plt.boxplot(highs, notch = False, showmeans = True, meanline = True)#, showfliers = False)
# plt.xticks([1, 2, 3], highnames)
# plt.ylabel('Energía libre (KT)')
# plt.show()

# plt.boxplot(lows, notch = False, showmeans = True, meanline = True)#, showfliers = False)
# plt.xticks([1, 2], lownames)
# plt.ylabel('Energía libre (KT)')
# plt.show()


# fig, axs = plt.subplots(1,5, sharey=True, sharex = True, tight_layout=True)


# axs[0].hist(NCeigvals[1], density = True, label = names[0])
# axs[1].hist(SL2eigvals[0], density = True, label = names[1])
# axs[2].hist(SL3eigvals[0], density = True, label = names[2])
# axs[3].hist(NCSL2eigvals[0], density = True, label = names[3])
# axs[4].hist(NCSL3eigvals[0], density = True, label = names[4])
# plt.ylim( 0, 0.6 )
# plt.legend()
# plt.show()

# E_NCSL3 = np.mean(NCSL3energy) - np.mean(NCenergy) - np.mean(SL3energy)
# E_NCSL2 = np.mean(NCSL2energy) - np.mean(NCenergy) - np.mean(SL2energy)





