%%
clc;clear;
addpath(genpath(pwd))

expirement_count = 52;  % change experiment number

config_sim;             % change paths and sim parameters
init_sim;               % Init sim - don't touch 
run_sim;                % main sim loop - don't touch 
thefigure_v4
% thefigure_v6