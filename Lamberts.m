function [V1, V2] = Lamberts(R1, R2, t)

muearth = 398600; % km^3/s^2

%% Prep
r1 = norm(R1);
r2 = norm(R2);
cross12   = cross(R1, R2);
Dtheta = acos(dot(R1,R2)/(r1*r2));

if cross12(3) < 0
    Dtheta = 2*pi - Dtheta;
end

A = sin(Dtheta)*sqrt(r1*r2/(1 - cos(Dtheta)));

% Determine approximately where F changes sign
% Use that value of z as the starting value
z = -100;
func = (y(z)/stumpC(z))^1.5*stumpS(z) + A*sqrt(y(z)) - sqrt(muearth)*t;
while func < 0
    z = z + 0.1;
    func = (y(z)/stumpC(z))^1.5*stumpS(z) + A*sqrt(y(z)) - sqrt(muearth)*t;
end


%% Lamberts

% Set an error tolerance and a limit on the number of iterations
TOL = 1.e-8;
maxItr = 5000;

% Iterate until z is determined to within the error tolerance
ratio = 1;
itr = 0;

while (abs(ratio) > TOL) && (itr <= maxItr)
    
    f = (y(z)/stumpC(z))^1.5*stumpS(z) + A*sqrt(y(z)) - sqrt(muearth)*t;
    
    if z == 0
        df = (sqrt(2)/40)*(y(0)^1.5) + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
    else
        df = (y(z)/stumpC(z))^1.5*(1/2/z*(stumpC(z) - 3*stumpS(z)/2/stumpC(z)) ...
               + 3*stumpS(z)^2/4/stumpC(z)) + A/8*(3*stumpS(z)/stumpC(z)*sqrt(y(z)) ...
               + A*sqrt(stumpC(z)/y(z)));
    end
    
    ratio = f/df;
    z = z - ratio;
    itr = itr + 1;
end

if itr >= maxItr
    fprintf('\n\n **Number of iterations exceeded** \n\n')
end


f = 1 - y(z)/r1;
g = A*sqrt(y(z)/muearth);
gdot = 1 - y(z)/r2;

V1 = (1/g)*(R2 - f*R1);
V2 = (1/g)*(gdot*R2 - R1);

function y = y(z)
    y = r1 + r2 + A*(z*stumpS(z) - 1)/sqrt(stumpC(z));
end


end