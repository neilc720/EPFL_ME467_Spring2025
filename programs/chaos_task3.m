close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
t_study = 2500.0;              % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
T_max = 120.0;                 % maximum T for recurrent flow analysis
T_eqb = 10.0;                  % time interval for computing equilibria
epsilon = 1e-6;                % perturbation amplitude for computing J

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% time integrate the KSE
[V,t_vec] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm);

%% computing equilibria
snapshots = [0,10,200,666,1111,2000];     % Modify the set of snapshot numbers. Do NOT
                               % use the rand or randi command to select
                               % snapshots from the trajectory; choose six 
                               % arbitrary numbers yourself. 

equilibria = {};                        
E_list = []; P_list = []; D_list = []; 
success_idx = [];                      

for i = 1:length(snapshots)
    idx = snapshots(i) + 1;            
    v_guess = V(:, idx);

    % Recherche de l'équilibre à partir de ce snapshot
    [v_eq, flag] = search4EQ(v_guess, T_eqb, dt, L, N, symm);

    if flag == 1
        fprintf(' Converged equilibrium at snapshot #%d (t = %.1f)\n', snapshots(i), t_vec(idx));
        success_idx(end+1) = snapshots(i);

        % Convertir en état physique
        u_guess = vector2field(v_guess, N, symm);
        u_eq = vector2field(v_eq, N, symm);

        % Tracer comparaison guess vs solution (Plot C)
        figure;
        plot(x, u_guess, 'k--', 'LineWidth', 1.5); hold on;
        plot(x, u_eq, 'r-', 'LineWidth', 2);
        xlabel('x'); ylabel('u(x)');
        legend('Initial guess','Equilibrium');
        title(sprintf('Plot C: Equilibrium from snapshot #%d', snapshots(i)));
        grid on;

        % Stockage
        equilibria{end+1} = v_eq;
        [E, P, D] = projection(v_eq, N, L, symm);
        E_list(end+1) = E;
        P_list(end+1) = P;
        D_list(end+1) = D;

    else
        fprintf(' No convergence at snapshot #%d (t = %.1f)\n', snapshots(i), t_vec(idx));
    end
end

[E_traj, P_traj, D_traj] = projection(V, N, L, symm);

figure;
plot3(E_traj, P_traj, D_traj, '-', 'Color', [0 0 0 0.2]); hold on;


scatter3(E_list, P_list, D_list, 60, 'r', 'filled');
xlabel('E'); ylabel('P'); zlabel('D');
title('Plot E: Attractor and equilibria');
grid on; view(3);
legend('Chaotic trajectory', 'Equilibria');

%% computing UPOs
% 4.2 - Compute recurrence indicator r(t,T)
guesses=[420 95; 330 115; 2191 85; 1603 103; 1066 57; 793 51];
T_vals = 5:1:T_max;  % valeurs de T à tester
nT = length(T_vals);
nt = size(V, 2);

R = zeros(nt, nT);     % matrice des distances de récurrence

for i = 1:nt
    vi = V(:,i);
    norm_vi = norm(vi);

    for j = 1:nT
        Tshift = round(T_vals(j)/dt_store);  % décalage en nombre de snapshots
        if (i + Tshift) <= nt
            vj = V(:, i + Tshift);
            R(i,j) = norm(vj - vi) / norm_vi;
        else
            R(i,j) = NaN;
        end
    end
end

%% Affichage de ln(r)
logR = log(R');
figure;
h = imagesc(T_vals, t_vec, logR);
set(h, 'AlphaData', ~isnan(logR));  % Make NaNs transparent

set(gca, 'YDir', 'normal');
xlabel('T'); ylabel('t');
colorbar;
title('Plot F: log(r(t,T)) — recurrence indicator');

upos = {}; T_list = []; lambda1_list = []; tF_list = [];
%%
for k = 1:size(guesses,1)
    t_guess = guesses(k,1);
    T_guess = guesses(k,2);

    idx = round(t_guess / dt_store) + 1;
    v_guess = V(:,idx);

    [v_po, T_po, flag] = search4PO(v_guess, T_guess, dt, L, N, symm);

    if flag == 1
        fprintf(' UPO converged at t = %.1f, T ≈ %.2f\n', t_guess, T_po);

        % Projection dans (E,P,D)
        [E,P,D] = projection(v_po, N, L, symm);
        upos{end+1} = v_po;
        T_list(end+1) = T_po;

        % Calcul du Floquet multiplier
        J = Jacobian(v_po, T_po, epsilon, dt, L, N, symm);
        lambda = eig(J);
        [lambda1, ~] = max(abs(lambda));
        lambda1_list(end+1) = lambda1;

        tF = T_po / log(lambda1);
        tF_list(end+1) = tF;

    else
        fprintf(' UPO did NOT converge at t = %.1f, T = %.2f\n', t_guess, T_guess);
    end
end


% Attracteur
[E_traj, P_traj, D_traj] = projection(V, N, L, symm);

figure;
plot3(E_traj, P_traj, D_traj, '-', 'Color', [0 0 0 0.2]); hold on;

% Ajouter les UPOs
for i = 1:length(upos)
    [E, P, D] = projection(upos{i}, N, L, symm);
    scatter3(E, P, D, 60, 'b', 'filled');
end

xlabel('E'); ylabel('P'); zlabel('D');
title('Plot G: Attractor and UPOs');
grid on; view(3);
legend('Trajectory','UPOs');
fprintf('\nUPO Summary Table:\n');
fprintf('%5s | %8s | %10s | %10s\n','#','T','|lambda1|','t*_F');
for i = 1:length(T_list)
    fprintf('%5d | %8.2f | %10.4f | %10.4f\n', i, T_list(i), lambda1_list(i), tF_list(i));
end
%% Finalizing 4.2: UPO integration, projection, stability, Plot G

% Prepare colormap for distinct colors
colors = lines(length(upos));

% Time-integrate each UPO over one full period to visualize the orbit shape
for i = 1:length(upos)
    vpo = upos{i};
    Tpo = T_list(i);

    % Integrate UPO over one period, storing ~100 snapshots
    [V_po, t_po] = KSE_integrate(vpo, Tpo, dt, Tpo / 100, L, N, symm);

    % Project trajectory in (E,P,D)
    [E_po, P_po, D_po] = projection(V_po, N, L, symm);

    % Overlay the orbit loop onto the attractor
    plot3(E_po, P_po, D_po, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    hold on;

    % Also mark starting point
    scatter3(E_po(1), P_po(1), D_po(1), 60, colors(i,:), 'filled');
end

% Final formatting of Plot G
xlabel('E'); ylabel('P'); zlabel('D');
title('Plot G: Attractor and UPOs with their orbits');
grid on; view(3);

% Build legend entries
legend_entries = ["Chaotic Trajectory"];
for i = 1:length(upos)
    legend_entries(end+1) = sprintf("UPO #%d", i);
end
legend(legend_entries, 'Location', 'best');

% Print summary table with Floquet multipliers and divergence times
fprintf('\nUPO Summary Table:\n');
fprintf('%5s | %8s | %10s | %10s\n', '#', 'T', '|lambda1|', 't*_F');
for i = 1:length(T_list)
    fprintf('%5d | %8.2f | %10.4f | %10.4f\n', ...
        i, T_list(i), lambda1_list(i), tF_list(i));
end
