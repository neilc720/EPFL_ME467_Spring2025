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
r = zeros(size(v1000));  
for k = 1:length(r)
    r(k) = epsilon * (2 * rand() - 1) * v1000(k);  
end

v0_perturbed = v1000 + r; 

[vv1, tt] = KSE_integrate(v1000, T_study, dt, dt_store, L, N, symm);
[vv2, ~]  = KSE_integrate(v0_perturbed, T_study, dt, dt_store, L, N, symm);

u1 = zeros(N, length(tt));
u2 = zeros(N, length(tt));
for j = 1:length(tt)
    u1(:,j) = vector2field(vv1(:,j), N, symm);
    u2(:,j) = vector2field(vv2(:,j), N, symm);
end

% Plot results
figure('Position', [100 100 1000 700])

subplot(3,1,1)
imagesc(tt, x, u1)
axis xy
xlabel('Time t'); ylabel('x')
title('Reference trajectory: u_1(x,t)')
colorbar

subplot(3,1,2)
imagesc(tt, x, u2)
axis xy
xlabel('Time t'); ylabel('x')
title('Perturbed trajectory: u_2(x,t)')
colorbar

subplot(3,1,3)
imagesc(tt, x, abs(u1 - u2))
axis xy
xlabel('Time t'); ylabel('x')
title('Difference |u_1(x,t) - u_2(x,t)|')
colorbar

sgtitle('Plot A: Sensitivity to Perturbation in KSE System')

