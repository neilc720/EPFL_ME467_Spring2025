close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 22.0;                      % domain length
N = 64;                        % spatial resolution
symm = false;                  % imposed center symmetry
dt = 0.1;                      % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2;                       % time intervals of computing J
N_norm = 5000;                 % number of re-normalizations
N_exp = 10;                    % number of Lyapunov exponents

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L) ...
    + cos(4.0*pi*x/L);         % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% to be completed...