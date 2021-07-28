function [dv,dt,ecc,h,a,T,rp,ra] = Hohmann(Orbit1, Orbit2)

% ===================================================
%
% Simple Inner to Outter Hohmann Transfer
%
%   Orbit1  - Inner Orbit
%   Orbit2  - Outter Orbit
%   dv      - Velocity needed for transfer
%   dt      - time it takes to transfer
%   ecc     - transfer eccentricity
%   h       - transfer angular momentum
%   a       - transfer semi-major axis
%   T       - transfer period
%
% ===================================================

MU  = 398600; % km^3/s^2

rp  = Orbit1.rp;
ra  = Orbit2.ra;
ecc = (ra - rp)/(ra + rp);

hi  = Orbit1.h;
h  = sqrt(rp*MU*(1 + ecc));
hf  = Orbit2.h;

Vi  = hi/rp;
Vti = h/rp;
Vtf = h/ra;
Vf  = hf/ra;
dv  = abs(Vi - Vti) + abs(Vtf - Vf);

a   = (ra + rp)/2;
T   = (2*pi)/sqrt(MU)*a^(3/2);
dt  = T/2;

end

