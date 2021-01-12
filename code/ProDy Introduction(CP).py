#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 22:12:37 2020
@author: aarondomenzain
"""


from prody import *
import matplotlib.pyplot as plt
import numpy as np
import os

#Compute free energy, eigenvalues and eigenvectors of a Gaussian Network Model from given PDB structure
    # Arguments:
        # protPDB; PDB ID of structure to apply GNM
        # model_number: Structure index (within the same PDB) to apply GNM
    # Outputs:
        # 3 element list: Free Energy (in KbT) , Eigenvalues, Eigenvectors

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
    
    # Gaussian Network Modelling
    gnm = GNM(protPDB)
    gnm.buildKirchhoff(network, cutoff = 10.0 , gamma = 1.0 )

    # Kirchoff Matrix
    K = gnm.getKirchhoff()

    # Calculate Normal Modes
    gnm.calcModes(1000) #All normal modes (default is 20)

    # Get Eigenvalues and eigenvectors
    Eigvals = gnm.getEigvals() #This is a vector containing eigenfrequencies
    Eigvecs = gnm.getEigvecs() #This is a matrix in which each column is an eigenvector (normal mode)
    

    # Generate mobility plot (Squared fluctuations in arbitrary units) for each mode for the FIRST model
    if model_number == 1:
        
        # First 5 modes
        n_modes = 5
    
        # plotmobility(protPDB, gnm, n_modes) 

    
    # Compute Free Energy F in KbT units 
    Kb = 1.380649 * (10 ** -23) # Boltzmann constant (in J/K)
    hbar = 3.16152649 * (10 ** -26) # Planck's reduced constant (in J*s)
    T = 310 # Human body temperature (in K)
    N = network.numAtoms() # Number of nodes in the network
    
    # We force the number of nodes to be equal to the number of modes for NC+SL3 complex
    # if protPDB is '1a1t':
    #     N = 74
    
    # Frecuencies w are proportional to the squared root of eigenvalues K
    frequencies = np.sqrt( Eigvals ) 
    
    E =  - 3 * (N-1) * np.log( (Kb * T) / hbar )  +  (1/1) * np.sum( np.log( frequencies ) )
    
    return (E, Eigvals, Eigvecs)


#Iterates over the PDB structure number model 
def gnmanalysis_models(protPDB, n_models):
    #Initialize a list of zeros of length (n_models)
    listEnergies = np.zeros(n_models)
    
    for n in range (n_models):
        
        Energy, Eigvals, Eigvecs = gnmanalysis(protPDB, n+1)
        
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
        
        # Checks if directory exists. If not, creates a directory with the PDB ID name
        if os.path.isdir('../results/' + protPDB + '/') is False:
            os.mkdir('../results/' + protPDB + '/')

    np.save('../results/' + protPDB + "/" + protPDB + "_model_energies.npy", listEnergies)
    np.save('../results/' + protPDB + "/" + protPDB + "_eigenvalues.npy", listEigvals)
    np.save('../results/' + protPDB + "/" + protPDB + "_eigenvectors.npy", listEigvecs)
    # # return (listEnergies, listEigvals, listEigvecs)




    #Iterates over a list of specific single-model PDB files associated with one biological structure
    
    # protPDB_list is a list containing custom pdb file names as strings in ' ' or " ". Termination .pdb must be specified
    # protPDB_name is a string containing the name of the protein asociated with all models
def gnmanalysis_custom(protPDB_list, protPDB_name):
    
    n = 0
    for protPDB in protPDB_list:
        Energy, Eigvals, Eigvecs = gnmanalysis(protPDB, 1)
        
        # No encontré manera bonita de iniciar matrices de ceros porque necesito conocer el número de eigenvectores y su tamaño
        if n == 0:
            
            N_models = len(protPDB_list)
            Neigvals = len(Eigvals) 
            Neigvecs = len(Eigvecs)
            
            #Calculate dimensions of eigenvectors matrix
            [N_vecs, vecs_length] = Eigvecs.shape
            
            #Initialize zeros array
            listEigvals = np.zeros([N_models, Neigvals])
            listEigvecs = np.zeros([N_models, N_vecs, vecs_length])
            listEnergies = np.zeros( N_models )
        
        listEnergies[n] = Energy
        listEigvals[n] = Eigvals
        listEigvecs[n] = Eigvecs
        
        n = n + 1
        # Checks if directory exists. If not, creates a directory with the PDB ID name
        if os.path.isdir('../results/' + protPDB_name + '/') is False:
            os.mkdir('../results/' + protPDB_name + '/')

    np.save('../results/' + protPDB_name + "/" + protPDB_name + "_model_energies.npy", listEnergies)
    np.save('../results/' + protPDB_name + "/" + protPDB_name + "_eigenvalues.npy", listEigvals)
    np.save('../results/' + protPDB_name + "/" + protPDB_name + "_eigenvectors.npy", listEigvecs)

# Draw mobility plot for first N modes
    # Arguments:
        # protPDB:PDB ID of desired protein
        # gnm: Gaussian Network Model ProDy instance
        # N: Number of modes to process
    
    # Outputs:
        # N images in .png format with name: protPDB/protPDB_SqFlucts_Mode_i.png with i= 1,2,..,N
    
def plotmobility(protPDB, gnm, N):
    # Checks if directory exists. If not, creates a directory with the PDB ID name
    if os.path.isdir('../results/' + protPDB + '/') is False:
        os.mkdir('../results/' + protPDB + '/')

    # First 5 lowest frequency modes 
    for i in range (N):
        if i < N:
            
            plt.figure()
            showSqFlucts(gnm[i], hinges=True)
            
            if i < 2:
                plt.ylim(0, 3)
            else:
                plt.ylim(0, 0.5)
                
            plt.ylabel("Square fluctuations (A. U.)")
            plt.xlabel("Residue")
            plt.savefig('../results/' + protPDB + '/' + protPDB + "_SqFlucts_Mode_" + str(i+1) + ".png", dpi=300)
            plt.show()
            
            
    
#####################      INSTRUCTIONS       ##########################
    
# GNM Normal Mode Analysis
    

# PDB list to analyze:
    #Unbounded 
#     '1aaf' : Unbounded NC from HIV-1 (20 models) color = 'gray'
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
    
# gnmanalysis_model

# gnmanalysis_models('1aaf', 20)
# gnmanalysis_models('1esy', 19)
# gnmanalysis_models('1bn0', 11)
# gnmanalysis_models('1f6u', 20)
# gnmanalysis_models('1a1t', 25)
# gnmanalysis_models('1a1t.pdb', 25)


# 4. LOAD DATA

# NCenergy = np.load("../results/1aaf/1aaf_model_energies.npy")
# SL2energy = np.load("../results/1esy/1esy_model_energies.npy")
# SL3energy = np.load("../results/1bn0/1bn0_model_energies.npy")
# NCSL2energy = np.load("../results/1f6u/1f6u_model_energies.npy")
# # NCSL3energy = np.load("../results/1a1t/1a1t_model_energies.npy")
# NCSL3energy = np.load("../results/1a1t.pdb/1a1t.pdb_model_energies.npy")

# NCeigvals = np.load("../results/1aaf/1aaf_eigenvalues.npy")
# SL2eigvals = np.load("../results/1esy/1esy_eigenvalues.npy")
# SL3eigvals = np.load("../results/1bn0/1bn0_eigenvalues.npy")
# NCSL2eigvals = np.load("../results/1f6u/1f6u_eigenvalues.npy")
# # NCSL3eigvals = np.load("../results/1a1t/1a1t_eigenvalues.npy")
# NCSL3eigvals = np.load("../results/1a1t.pdb/1a1t.pdb_eigenvalues.npy")

# NCsim_energy = np.load( "../results/NC_simulation/NC_simulation_model_energies.npy" ) 
# NCsim_eigvals = np.load( "../results/NC_simulation/NC_simulation_eigenvalues.npy" ) 

# NCSL2sim_energy = np.load( "../results/NC-SL2_simulation/NC-SL2_simulation_model_energies.npy" ) 
# NCSL2sim_eigvals = np.load( "../results/NC-SL2_simulation/NC-SL2_simulation_eigenvalues.npy" ) 


# # Calculate binding energy from mean free energy changes
# E_NCSL3 = np.mean(NCSL3energy) - np.mean(NCenergy) - np.mean(SL3energy)
# E_NCSL3_std = ( np.std(NCSL3energy) - np.std(NCenergy) - np.std(SL3energy) ) / np.sqrt( NCSL3eigvals.shape[1] )


# E_NCSL2 = np.mean(NCSL2energy) - np.mean(NCenergy) - np.mean(SL2energy)
# E_NCSL2_std = ( np.std(NCSL2energy) - np.std(NCenergy) - np.std(SL2energy) ) / np.sqrt( NCSL2eigvals.shape[1] )

# # ztest = ( np.mean( NCSL3energy ) - np.mean( NCSL2energy ) ) / np.sqrt(( np.std(NCSL3energy)**2 - np.std(NCSL2energy)**2 ))

# print()
# print("NC-SL2 binding energy :", E_NCSL2)
# print("NC-SL3 binding energy :", E_NCSL3)


# # ## Free energy distributions
# #     # Labels for plots
# complex = [SL2energy, SL3energy, NCenergy, NCSL2energy, NCSL3energy]
# names = ["SL2", "SL3", 'NC', 'NC+SL2', 'NC+SL3']

# highs = [NCenergy, NCSL2energy, NCSL3energy]
# highnames = ['NC', 'NC+SL2', 'NC+SL3']

# SLcomplexes = [NCSL2energy, NCSL3energy]
# SLcomplexes_names = ['NC+SL2', 'NC+SL3']

# plt.boxplot(complex, notch = False, showmeans = True, meanline = True, showfliers = False)
# plt.xticks([1, 2, 3, 4, 5], names)
# plt.ylabel('Energía libre (KT)')
# plt.show()
    
# plt.boxplot(highs, notch = False, showmeans = True, meanline = True, showfliers = False)
# plt.xticks([1, 2, 3], highnames)
# plt.ylabel('Energía libre (KT)')
# plt.show()

# plt.boxplot(SLcomplexes, notch = False, showmeans = True, meanline = True, showfliers = True)
# plt.xticks([1, 2], SLcomplexes_names)
# plt.ylabel('Energía libre (KT)')
# plt.show()







# fig, axs = plt.subplots(1,3, sharey=False, sharex = True, tight_layout=True)


# axs[0].hist(NCeigvals[1], density = True, label = names[0])
# axs[1].hist(BCSL2eigvals[0], density = True, label = names[1])
# axs[2].hist(NCSL3eigvals[0], density = True, label = names[2])
# plt.ylim( 0, 0.6 )
# plt.legend()
# plt.show()








