function coes = COES(R,V)

% =================================================
%
% COES : Obtains classical orbital elements
%
%   R - Position Vector
%   V - Velocity Vector
%
%   h     - angular momentum
%   i     - inclination
%   RAAN  - right ascension of ascending node
%   e     - eccentricity
%   w     - argument of perigee
%   TA    - true anomaly
%   T     - period
%   t     - current time on period
%   a     - semi-major axis
%   rp    - radius of perigee
%   ra    - radius of apogee
%   Me    - Mean Anomaly
%   EE    - Eccentric Anomaly
%
% =================================================


TOL  = 1.e-10;
MU   = 398600;         % km^3/s^2

r    = norm(R);      
v    = norm(V);        
vr   = dot(R, V)/r;

H    = cross(R, V);  
h    = norm(H);   

inc  = acos(H(3)/h);  

N    = cross([0 0 1], H);
n    = norm(N);

RAAN = raan(N, n);

E    = (1/MU)*(cross(V, H) - MU*(R/r));
ecc  = norm(E);

w    = omega(E, ecc, N, n, TOL);

TA   = theta(E, ecc, N, n, R, r, vr, TOL);

a    = h^2/MU/(1 - ecc^2);
ra   = h^2/MU/(1 - ecc);
rp   = h^2/MU/(1 - ecc);

EE   = ee(ecc, TA);
Me   = EE - ecc*sin(EE);

T    = 2*pi*a^(3/2)/sqrt(MU);
t    = abs(Me*T/(2*pi));


coes = [h, inc, RAAN, ecc, w, TA, T, t, a, ra, rp, Me, EE];

function temp = raan(N, n)
    if n ~= 0
        temp = acos(N(1)/n);
        if N(2) < 0
            temp = 2*pi - temp;
        end
    else
        temp = 0;
    end
end

function temp = omega(E, ecc, N, n, TOL)
    if n ~= 0
        if ecc > TOL
           temp = acos(dot(N,E)/(n*ecc));
            if E(3) < 0
                temp = 2*pi - temp;
            end
        else
            temp = 0;
        end
    else
        temp = 0;
    end
end

function temp = theta(E, ecc, N, n, R, r, vr, TOL)
    if ecc > TOL
        temp = acos(dot(E,R)/(ecc*r));
        if vr < 0
            temp = 2*pi - temp;
        end
    else
        cp = cross(N,R);
        if cp(3) >= 0
            temp = acos(dot(N,R)/(n*r));
        else
            temp = 2*pi - acos(dot(N,R)/(n*r));
        end
    end
end

function temp = ee(ecc, TA)
    if TA > pi
        temp = atan(sqrt((1 - ecc)/(1 + ecc))*tan((TA-pi)/2))*2 + pi;
    else
        temp = atan(sqrt((1 - ecc)/(1 + ecc))*tan(TA/2))*2;
    end
end


end

