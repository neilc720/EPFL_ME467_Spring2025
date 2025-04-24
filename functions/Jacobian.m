%%% DESCRIPTION -----------------------------------------------------------
%   Jacobian of the KSE flow map: J=dv(t)/dv(0)


%%% INPUTS ----------------------------------------------------------------
%   v0          reference point of the Jacobian (column state vector)
%   t           integration time interval
%   epsilon     perturbation magnitude for finite difference derivatrives
%   dt          step size in time integrations
%   L           domain length
%   N           spatial resolution
%   symm        center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   J       Jacobian matrix


function J = Jacobian(v0,t,epsilon,dt,L,N,symm)
    n = length(v0);
    J = zeros(n,n);
    
    %%% to be completed...
end
