% a simplified steady state model for populatin density change of pairwise coculture
function [df,N1ss,N2ss] = mse_2mer_ss(log10z,N0,G0,FC)

% Input variables:
%       log10z:     log10(multiplication of Yield_aa, Rcarbon and Byp_frac)
%       N0:         initial cell density
%       G0:         initial glucose concentration

% Output variables:
%       df:         difference between observed and simulated fold changes
%       N1, N2:     steady state cell density

z= 10.^log10z;

if (z(1)>z(2))
    delta1 = N0(1) + N0(2) + G0 * z(2);
    delta2 = z(1) - z(2);
    delta3 = (G0 * z(2) + 2 * N0(2)) * G0;
    N1ss     = N0(1) - (delta1 - sqrt(delta1^2 + delta2*delta3)) / delta2 * z(1);
    N2ss     = N0(2) + (G0 + (N0(1)-N1ss)/z(1)) * z(2);
else
    delta1 = N0(1) + N0(2) + G0 * z(1) ;
    delta2 = z(2) - z(1);
    delta3 = (G0 * z(1) + 2 * N0(1)) * G0 ;
    N2ss     = N0(2) - (delta1 - sqrt(delta1^2 + delta2*delta3)) / delta2 * z(2);
    N1ss     = N0(1) + (G0 + (N0(2)-N2ss)/z(2)) * z(1);
end

df = [N1ss/N0(1)-FC(1); N2ss/N0(2)-FC(2)];

end

