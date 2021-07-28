function approx = newtonsMethod(f, fprime, x0, TOL, Me, ecc)

% ==========================================
%
% Modified version of Newton's Iterative Method
%
%   f      - Me - E + ecc*sin(E);
%   fprime - -1 + ecc*cos(E);
%   x0     - Inital guess
%   TOL    - Tolerance
%   Me     - Mean Anomaly
%   ecc    - eccentricity
%
% ==========================================


MAX_ITR = 12;
x2      = x0 - f(Me, x0, ecc)/fprime(x0, ecc);
itr     = 1;
error   = abs( x2 - x0 );

while error > TOL && itr < MAX_ITR
    x1    = x2;
    x2    = x1 - f(Me, x1, ecc)/fprime(x1, ecc);
    error = abs(x2 - x1);
    itr   = itr + 1;   
end

approx = x2;
end