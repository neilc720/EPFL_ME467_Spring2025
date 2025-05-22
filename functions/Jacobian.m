function J = Jacobian(v0, t, epsilon, dt, L, N, symm)
    % v0      : initial state vector
    % t       : integration time (scalar)
    % epsilon : perturbation magnitude
    % dt      : integration time step
    % L, N    : domain parameters
    % symm    : symmetry flag (boolean)

    n = length(v0);        % size of the state vector
    J = zeros(n, n);       % allocate Jacobian matrix

    % Loop over each component of v0
    for i = 1:n
        % Unit vector in direction i
        ei = zeros(n, 1);
        ei(i) = 1;

       
        v_plus  = v0 + epsilon * ei;
        v_minus = v0 - epsilon * ei;

       
        [v_final_plus, ~]  = KSE_integrate(v_plus,  t, dt, 0, L, N, symm);
        [v_final_minus, ~] = KSE_integrate(v_minus, t, dt, 0, L, N, symm);

        % Central difference approximation for column i
        J(:, i) = (v_final_plus - v_final_minus) / (2 * epsilon);
    end
end
