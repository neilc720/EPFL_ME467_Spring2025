close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
T_trans = 1000.0;              % transient time period
T_study = 250.0;               % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
epsilon = 1e-2;                % relative ampliture of perturbation

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% transient time integration
[v1000,~] = KSE_integrate(v0,T_trans,dt,0,L,N,symm);

%% perturbing the state vector
r = zeros(size(v1000));        % the perturbation vector memory allocation
for k = 1:length(r)
    %%% to be completed
end

%%% to be completed