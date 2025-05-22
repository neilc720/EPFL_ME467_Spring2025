close all; clear; clc

%% Inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
dt = 0.1;                      % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2;                       % time interval between Jacobian applications
N_norm = 5000;                 % number of re-normalizations
N_exp = 10;                    % number of Lyapunov exponents to compute

%% Initial condition (Eq. 17)
[x,~] = domain(L,N);           
u0 = sin(2.0*pi*x/L);          
v0 = field2vector(u0,N,symm);  

%% Initialize orthonormal deviation vectors (W₀ = identity columns)
Q = zeros(length(v0), N_exp);
for i = 1:N_exp
    Q(i,i) = 1;
end

%% Allocate memory
log_sum = zeros(N_exp, 1);     % Cumulative log(|Rjj|) for each exponent
X = zeros(N_exp, N_norm);      % Running average of exponents
t = zeros(N_norm, 1);          % Time points

%% Main loop
figure
for i = 1:N_norm
    % Step 1: Jacobian of the flow map from t=0 to t=tau
    J = Jacobian(v0, tau, epsilon, dt, L, N, symm);

    % Step 2: Apply Jacobian to current orthonormal vectors
    V = J * Q;

    % Step 3: QR decomposition (V = Q R)
    [Q, R] = qr(V, 0);

    % Step 4: Accumulate log-growth from R diagonal
    for j = 1:N_exp
        log_sum(j) = log_sum(j) + log(abs(R(j,j)));    % ← X ← X + ln|Rjj|
        X(j, i) = log_sum(j) / (i * tau);              % χ_j(t) = average
    end

    % Step 5: Advance reference trajectory
    [v0, ~] = KSE_integrate(v0, tau, dt, 0, L, N, symm);
    t(i) = i * tau;

    % Step 6: Plot every 20 renormalizations
    if mod(i, 20) == 0
        clf; hold on; grid on
        colors = lines(N_exp);
        legend_entries = cell(N_exp, 1);
        for q = 1:N_exp
            plot(t(1:i), X(q,1:i), 'LineWidth', 1.5, 'Color', colors(q,:));
            legend_entries{q} = ['$\chi_', num2str(q), '$'];
        end
        xlabel('Time $t$', 'Interpreter', 'latex')
        ylabel('$\chi_i(t)$', 'Interpreter', 'latex')
        title('Plot B: Convergence of Lyapunov Exponents', 'Interpreter', 'latex')
        legend(legend_entries, 'Interpreter', 'latex', 'Location', 'eastoutside')
        drawnow
    end
end

%% Final values
chi_final = X(:, end);
disp('Final Lyapunov exponents (χ₁ to χ₁₀):');
disp(chi_final);
lyapunov_time = 1/max(chi_final)
