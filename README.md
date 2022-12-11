# rl-qlearn
Supplementary material for "A residual life value function for dynamic treatment regime optimization via Q-learning"

This repository contains the R and Python code used to conduct the simulation study in "A residual life value function for dynamic treatment regime optimization via Q-learning" by Rhodes, Davidian, and Lu.

To facilitate running the provided programs, we recommend creating the directory structure depicted in 'directory_structure.png' and saving all provided code in the parent directory 'Research.'

The program files should be executed in numerical order according to their file names. 

1) 01_generate_simulation_data.R - Generates the data sets for the simulation study.
2) 02_create_simulation_contextVecs.py - Creates the context vectors for the simulation data, and saves the trained LSTM autoencoders.
3) 03_generate_baseline_validation_data.R - Generates baseline data for the validation procedure.
4) 04_create_validation_contextVecs.py - Contains functions to create the context vectors for the validation data.
5) 05_conduct_qLearning.R - Conducts RL Q-learning to estimate the testing value and validation value of a user-specified treatment regime. 
6) 06_conduct_validation_obsNoTrt.R - Computes the validation value estimate for the observed treatment regime and the no treatment regime.
7) functions.R - Contains functions used to conduct the RL Q-learning method. 
8) sim_functions.R - Contains functions used to generate the data for the simulation study. 
