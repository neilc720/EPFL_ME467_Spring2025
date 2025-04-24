%%% DESCRIPTION -----------------------------------------------------------
%   dealiasing a 1D periodic field in spectral state using the 2/3 rule


%%% INPUTS ----------------------------------------------------------------
%   u   vector of Fourier coefficients (vector of complex numbers)


%%% OUTPUTS ---------------------------------------------------------------
%   u_  the same as 'u' with 1/3 of the highest frequency modes set to 0


function u_ = dealiase(u)
    N = length(u);
    Nd = round(N/3);
    
    u_ = u;
    u_(Nd+2:N-Nd) = 0;
end
