
% AERO 351
% PROJECT 1
% Isaac Velasco

clc ; close all;


%% Project Specifications

%   **Orbital Clean Up Mission
%
%       - Encounter Four Debris Objects
%       - Deliver Propulsion Package to each Object
%       - Must Stay with Object for 5 periods
%         to ensure Package has been delivered
%
%   **Must Use Four Types of Transfers:
%
%       - inc and/or RAAN change
%       - Phasing
%       - Hohmann
%       - Non-Hohmann
%       - Lambert
%       - Apse-Line
%       - Bi-Elliptical
%       - Non-Impulsive
%
%   **Asumme no Perturbations
%   **Total Delta-V must be < 18 km/s
%
%   **Two Objects in LEO
%   **One Object in MEO
%   **One Object in GEO
%       - Must be close to 0 ecc
%       - 0 < inc < 5 degrees
%

%   **Note: may be easier to choose all Objects with
%           prograde orbits


%% Orbit Specifications

%   **Hohmann
%       - Find speed at apoapse/periapse of initial orbit
%       - Find ecc of transfer orbit
%       - Find h of transfer orbit
%       - Find speed at apogee and perigee of transfer orbit
%       - Find speed at periapse/apoapse of final orbit
%       - Find difference of speed in both locations
%       **Note: If NOT circular orbit, transfer must be done at
%               periapse or apoapse
%             - Trajectory and transfer must share the same apse-line
%
%   **Bi-Elliptical
%       - Use if Radius of outer orbit > 15 times inner orbit
%       - Same as Hohmann, but need one extra burn and two transfer
%         ellipses
%
%   **Phasing
%       - Two impulse
%       - Inner Transfer
%           *Chaser behind target
%           *Can return to same location in less than one period
%               **Note: not reccomended
%       - Outer Transfer
%           *Chaser ahead of target
%           *Can return to same location in less than one period
%               **Note: not reccomended
%
%   **Non-Hohman
%       - Burns not done at apogee or perigee
%       - Must have SAME Apse-Line
%       - Need Magnitude & Direction
%       - Given Ra, Rp, TA(a), TA(b)
%           *Find ecc
%           *Find h
%           *Find V1 & V2
%           *Use Law of Cosines
%
%   **Apse-Line Rotation
%       - Two intersecting orbits
%       - Share same focus
%       - Steps
%           *Set rI1 = rI2
%           *Move h only terms to left
%           *Move other to right
%           *Use trig
%           *Sub Eta = TA1 - TA2
%           *Combine cos and sine terms
%
%   **RAAN and/or INC change
%       - Need common focus
%       - Maneuver occurs at apogee
%       - Use alpha for angle at which velocity is rotated
%
%   **Lambert
%
%   **Non-Impulsive
%

%   **NOTE: All equations found in Lectures 11-15


%% Process

%   **Use TLE's to acquire COE's for specified Object
%   **Use Mean Anomaly to solve for initial Eccentric Anomaly
%   **Use Newtons to solve for Final Eccentric Anomaly
%   **Use Final Eccentric to solve for True Anomaly
%   **Use previous information to solve for Current Position on orbit
%   **Use previous information to solve for Current Velocity on orbit
%   **Use Perifocal Frame to switch to ECI Frame
%   **Use vector format to propage 5 orbital periods necessary
%         to deliver package
%   **Use vector format to propagate to optimal orbit transfer 
%         locations (aka periapse and apoapse)
%   **Choose optimal orbit transfer type
%       - test various types if necessary


%% Constants & Useful Equations

REARTH = 6378;   % km
MU     = 398600; % km^3/s^2
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% COES
T      = @(rev) 86400/rev;
A      = @(ra, rp) (ra + rp)/2;          % semi-major axis
A2     = @(T) ((T^2*MU)/(2*pi)^2)^(1/3); % semi-major axis
H      = @(e, rp) sqrt(rp*MU*(1 + e));   % angular momentum
ECC    = @(ra, rp) (ra - rp)/(ra + rp);  % eccentricity
RA     = @(a, rp) 2*a - rp;              % radius of apogee
THETA  = @(E, e) atan(tan(E/2)*sqrt((1+e)/(1-e)))*2; % true anomaly

% Perifocal
RX     = @(h, e, TA) (h^2/MU)*(1/(1+e*cos(TA))).*[cos(TA); sin(TA); 0]; % km
VX     = @(h, e, TA)  MU/h .* [-sin(TA); e + cos(TA); 0];                % km/s

% ECI
R3     = @(n) [cos(n) sin(n) 0;
              -sin(n) cos(n) 0;
               0 0 1];
R1     = @(n) [1 0 0;
               0  cos(n) sin(n);
               0 -sin(n) cos(n)];           
QXX    = @(omega, inc, RAAN) R3(omega)*R1(inc)*R3(RAAN);
R_ECI  = @(Qxx, rx) transpose(Qxx)*rx; % km
v_ECI  = @(Qxx, vx) transpose(Qxx)*vx; % km/s


% Newton's Iteration
TOL    = 1.e-8;
F      = @(Me, E, ecc) Me - E + ecc*sin(E);
FPRIME = @(E, ecc) -1 + ecc*cos(E);

% Alpha for RAAN Change
ALPHAR  = @(inc, D_RAAN) acos(cos(inc)^2 + (sin(inc)^2)*cos(D_RAAN));
% Alpha for Inc and RAAN Change
ALPHARI = @(inc1, inc2, D_RAAN) acos(cos(inc1)*cos(inc2) + sin(inc1)*sin(inc2)*cos(D_RAAN));
% RAAN change deltav
DVRAAN  = @(alpha, V) 2*V*sin(alpha/2);

% Inclination no speed change
INOV    = @(dha, V) 2*V*sin(dha/2);
% Speed and inclination change congruently
INV     = @(dha, V1, V2) sqrt(V1^2 + V2^2 - 2*V1*V2*cos(dha/2));
% Plane change then speed change
VTHENI  = @(dha, V1, V2) 2*V1*sin(dha/2) + abs(V2-V1);
% Speed change then plane change
ITHENV  = @(dha, V1, V2) 2*V2*sin(dha/2) + abs(V2-V1);

% Circularization at Perigee
CHELPP  = @(ra, rp) sqrt(2*(-MU/(ra + rp) + MU/rp));
CIRCP   = @(rp) sqrt(2*(-MU/(2*rp) + MU/rp));


%% Launch Date

% November 14, 2020


%% GOAL

DELTAV = 0; % km/s


%% Objects

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~LEO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ----------------------Orbit Parameters-------------------------------   
% deb Transtage 4
%                   Epoch: Year - day + fractional portion                
% 1  3494U 65082QQ  20321.60609027  .00000162  00000-0  25246-4 0  9991
%          INC       RAAN    ECC     OMEGA    M        N
% 2  3494  31.8925 291.6972 0028451  25.1705 335.0316 14.58519088907316
%----------------------------------------------------------------------

LEOF1.T     = T(14.58519088);        % s
LEOF1.i     = deg2rad(31.8925);      % rads
LEOF1.ra    = 718 + REARTH;          % km
LEOF1.rp    = 678 + REARTH;          % km
LEOF1.a     = A(LEOF1.ra, LEOF1.rp); % km
LEOF1.e     = 0.0028451;               
LEOF1.h     = H(LEOF1.e, LEOF1.rp);  % km^2/s
LEOF1.RAAN  = deg2rad(291.6972);     % rads
LEOF1.w     = deg2rad(25.1705);      % rads
LEOF1.Me    = deg2rad(335.0316);     % rads


% ----------------------Orbit Parameters-------------------------------
% deb Transtage 4 
%                   Epoch: Year - day + fractional portion                
% 1  1964U 65082GZ  20321.58527382  .00000297  00000-0  10032-3 0  9996
%          INC       RAAN     ECC     OMEGA     M        N
% 2  1964  32.1892   289.9154 0193160 267.1062  90.7451  14.20152402827766
%----------------------------------------------------------------------

LEOF2.T     = T(14.20152402);        % s
LEOF2.i     = deg2rad(32.1892);      % rads
LEOF2.ra    = 964 + REARTH;          % km
LEOF2.rp    = 686 + REARTH;          % km
LEOF2.a     = A(LEOF2.ra, LEOF2.rp); % km
LEOF2.e     = 0.0193160;               
LEOF2.h     = H(LEOF2.e, LEOF2.rp);  % km^2/s
LEOF2.RAAN  = deg2rad(289.9154);     % rads
LEOF2.w     = deg2rad(267.1062);     % rads
LEOF2.Me    = deg2rad(90.7451);      % rads


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~MEO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ----------------------Orbit Parameters-------------------------------
%                                       T       INC     Za      Zp      ECC     
% 1979-101G  	29004	SATCOM 3 DEB	794.59	8.80	35,680	8,400	0.4799853
%                   Epoch: Year - day + fractional portion                      
% 1 29004U 79101G   20318.64683831 -.00000085  00000-0  00000-0 0  9990
%          INC     RAAN     ECC     OMEGA     M        N       
% 2 29004   8.7979 297.2829 4799853  63.8867 159.5878  1.81224489 97748
%----------------------------------------------------------------------

MEO.T     = 794.59*60;         % s
MEO.i     = deg2rad(8.80);     % rads
MEO.ra    = 35680 + REARTH;    % km
MEO.rp    = 8400 + REARTH;     % km
MEO.a     = A(MEO.ra, MEO.rp); % km
MEO.e     = 0.4799853;               
MEO.h     = H(MEO.e, MEO.rp);  % km^2/s
MEO.RAAN  = deg2rad(297.2829); % rads
MEO.w     = deg2rad(63.8867);  % rads
MEO.Me    = deg2rad(159.5878); % rads


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~GEO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ----------------------Orbit Parameters-------------------------------
%                                       T           INC     Za      Zp      ECC     
% 2018-050D  	43645	FENGYUN 2H DEB	1,436.36	0.31	35,987	35,596	0.0046335
%                   Epoch: Year - day + fractional portion                      
% 1 43645U 18050D   20316.59223193 -.00000264  00000+0  00000+0 0  9993
%           INC    RAAN     ECC      OMEGA    M        N  
% 2 43645   0.3138 265.4847 0046335  28.4479  65.6190  1.00253570  8424
%----------------------------------------------------------------------

GEO.T     = 1436.36*60;        % s
GEO.i     = deg2rad(0.31);     % rads
GEO.ra    = 35987 + REARTH;    % km
GEO.rp    = 35596 + REARTH;    % km
GEO.a     = A(GEO.ra, GEO.rp); % km
GEO.e     = 0.0046335;               
GEO.h     = H(GEO.e, GEO.rp);  % km^2/s
GEO.RAAN  = deg2rad(265.4847); % rads
GEO.w     = deg2rad(28.4479);  % rads
GEO.Me    = deg2rad(65.6190);  % rads



%% Initial State Vectors

% Initial Eccentric Anomaly Guess
E0L1 = EGuess(LEOF1.Me, LEOF1.e);
E0L2 = EGuess(LEOF2.Me, LEOF2.e);
E0M  = EGuess(MEO.Me, MEO.e);
E0G  = EGuess(GEO.Me, GEO.e);

% Final Eccentric Anomaly
EL1 = newtonsMethod(F, FPRIME, E0L1, TOL, LEOF1.Me, LEOF1.e);
EL2 = newtonsMethod(F, FPRIME, E0L2, TOL, LEOF2.Me, LEOF2.e);
EM  = newtonsMethod(F, FPRIME, E0M, TOL, MEO.Me, MEO.e);
EG  = newtonsMethod(F, FPRIME, E0G, TOL, GEO.Me, GEO.e);

% Get True Anomaly
LEOF1.TA = THETA(EL1, LEOF1.e); % rads
LEOF2.TA = THETA(EL2, LEOF2.e); % rads
MEO.TA   = THETA(EM, MEO.e);    % rads
GEO.TA   = THETA(EG, GEO.e);    % rads

% Perifocal Vectors
R0L1 = RX(LEOF1.h, LEOF1.e, LEOF1.TA); % km
R0L2 = RX(LEOF2.h, LEOF2.e, LEOF2.TA); % km
R0M  = RX(MEO.h, MEO.e, MEO.TA);       % km
R0G  = RX(GEO.h, GEO.e, GEO.TA);       % km

V0L1 = VX(LEOF1.h, LEOF1.e, LEOF1.TA); % km/s
V0L2 = VX(LEOF2.h, LEOF2.e, LEOF2.TA); % km/s
V0M  = VX(MEO.h, MEO.e, MEO.TA);       % km/s
V0G  = VX(GEO.h, GEO.e, GEO.TA);       % km/s

% ECI State Vectors
QL1 = QXX(LEOF1.w, LEOF1.i, LEOF1.RAAN);
QL2 = QXX(LEOF2.w, LEOF2.i, LEOF2.RAAN);
QM  = QXX(MEO.w, MEO.i, MEO.RAAN);
QG  = QXX(GEO.w, GEO.i, GEO.RAAN);

LEOF1.R = R_ECI(QL1, R0L1); % km
LEOF2.R = R_ECI(QL2, R0L2); % km
MEO.R   = R_ECI(QM, R0M);   % km
GEO.R   = R_ECI(QG, R0G);   % km
SC.R    = LEOF1.R;          % km

LEOF1.V = v_ECI(QL1, V0L1); % km/s
LEOF2.V = v_ECI(QL2, V0L2); % km/s
MEO.V   = v_ECI(QM, V0M);   % km/s
GEO.V   = v_ECI(QG, V0G);   % km/s
SC.V    = LEOF1.V;          % km/s


%% 1st LEO Object - 5 PERIODS

timespan = [0 LEOF1.T*5]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R', SC.V'];

% LEOF1 prop
[t1, S1] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [S1(end,1)  S1(end,2)  S1(end,3)];
LEOF2.R = [S1(end,7)  S1(end,8)  S1(end,9)];
MEO.R   = [S1(end,13) S1(end,14) S1(end,15)];
GEO.R   = [S1(end,19) S1(end,20) S1(end,21)];
SC.R    = LEOF1.R; 

LEOF1.V = [S1(end,4)  S1(end,5)  S1(end,6)];
LEOF2.V = [S1(end,10) S1(end,11) S1(end,12)];
MEO.V   = [S1(end,16) S1(end,17) S1(end,18)]; 
GEO.V   = [S1(end,22) S1(end,23) S1(end,24)];
SC.V    = LEOF1.V; 

% Check time based position along orbit, prop to apogee
C1L1 = COES(LEOF1.R,LEOF1.V);
timespan = [0 (LEOF1.T - C1L1(8) + LEOF1.T/2)];

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R', SC.V'];

% Prop to apogee for inclination change
[t1a, S1a] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [S1a(end,1)  S1a(end,2)  S1a(end,3)];
LEOF2.R = [S1a(end,7)  S1a(end,8)  S1a(end,9)];
MEO.R   = [S1a(end,13) S1a(end,14) S1a(end,15)];
GEO.R   = [S1a(end,19) S1a(end,20) S1a(end,21)];

LEOF1.V = [S1a(end,4)  S1a(end,5)  S1a(end,6)];
LEOF2.V = [S1a(end,10) S1a(end,11) S1a(end,12)];
MEO.V   = [S1a(end,16) S1a(end,17) S1a(end,18)]; 
GEO.V   = [S1a(end,22) S1a(end,23) S1a(end,24)];


%% LEO 1 Animation

% clear SA D vid
% SA = ones(length(S1(:,1))+length(S1a(:,1)),30);
% for i = 1:30
%     tmp = {[S1(:,i)],[S1a(:,i)]};
%     SA(:,i) = vertcat(tmp{:});
% end
% 
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% 
% set(gca,'XLim',[min([min(SA(:,1)) min(SA(:,7)) min(SA(:,13)) min(SA(:,19))])...
%                 max([max(SA(:,1)) max(SA(:,7)) max(SA(:,13)) max(SA(:,19))])],...
%         'YLim',[min([min(SA(:,2)) min(SA(:,8)) min(SA(:,14)) min(SA(:,20))])...
%                 max([max(SA(:,2)) max(SA(:,8)) max(SA(:,14)) max(SA(:,20))])],...
%         'ZLim',[min([min(SA(:,3)) min(SA(:,9)) min(SA(:,15)) min(SA(:,21))])...
%                 max([max(SA(:,3)) max(SA(:,9)) max(SA(:,15)) max(SA(:,21))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(307, 60);
% hold on
% 
% 
% for i = 1:6:length(SA(:,1))
%     addpoints(curve1,SA(i,1), SA(i,2), SA(i,3));
%     addpoints(curve2,SA(i,7), SA(i,8), SA(i,9));
%     addpoints(curve3,SA(i,13),SA(i,14),SA(i,15));
%     addpoints(curve4,SA(i,19),SA(i,20),SA(i,21));
%     
%     spccrft = scatter3(SA(i,1),  SA(i,2),  SA(i,3),'w>');
%     head2   = scatter3(SA(i,7),  SA(i,8),  SA(i,9), 'MarkerFaceColor', 'r','MarkerEdgeColor','r');
%     head3   = scatter3(SA(i,13), SA(i,14), SA(i,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c');
%     head4   = scatter3(SA(i,19), SA(i,20), SA(i,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m');
%     title('LEO 1')
%     legend('LEO_1','LEO_2','MEO','GEO','S/C & TRANSTAGE_4_ _D_E_B',...
%            'TRANSTAGE_4_ _D_E_B',...
%            'SATCOM_3_ _D_E_B','FENGYUN_2_H_ _D_E_B')
%     legend('TextColor','white')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(spccrft);
%     delete(head2);
%     delete(head3);
%     delete(head4);
% end
% 
% vid = VideoWriter('LEO1.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:6:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


% plot3(S2(:,1),  S2(:,2),  S2(:,3), 'b',...
%       S2(:,7),  S2(:,8),  S2(:,9), 'r',...
%       S2(:,13), S2(:,14), S2(:,15),'c',...
%       S2(:,19), S2(:,20), S2(:,21),'m')
% scatter3(S2(end,1),  S2(end,2),  S2(end,3), 'MarkerFaceColor', 'b','MarkerEdgeColor','b')
% scatter3(S2(end,7),  S2(end,8),  S2(end,9), 'MarkerFaceColor', 'r','MarkerEdgeColor','r')
% scatter3(S2(end,13), S2(end,14), S2(end,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c')
% scatter3(S2(end,19), S2(end,20), S2(end,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m')
% scatter3(S2(end,25), S2(end,26), S2(end,27),'k>')


%% 1st Transfer - inc Change/Circularize/Hohmann->Circularized LEO 2 ra

% Inclination Change at apogee
DELTAV = DELTAV + INOV(abs(LEOF1.i - LEOF2.i), norm(LEOF1.V));
SC     = LEOF1;  
SC.i   = LEOF2.i;

% Update COES before prop
TCTO   = COES(SC.R, SC.V);
SC.Me  = TCTO(12);

% Initial Eccentric Anomaly Guess
E0SC  = EGuess(SC.Me, SC.e);

% Final Eccentric Anomaly
ESC  = newtonsMethod(F, FPRIME, E0SC, TOL, SC.Me, SC.e);

% Get True Anomaly
SC.TA = THETA(ESC, SC.e); % rads

% Perifocal Vectors
R0SC   = RX(SC.h, SC.e, SC.TA); % km
V0SC   = VX(SC.h, SC.e, SC.TA); % km/s

% ECI State Vectors
Qsc    = QXX(SC.w, SC.i, SC.RAAN);
SC.R   = R_ECI(Qsc, R0SC); % km
SC.V   = v_ECI(Qsc, V0SC); % km/s

% Prop to perigee
timespan = [0 LEOF1.T/2];

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R, SC.V];
[tT1, T1] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [T1(end,1)  T1(end,2)  T1(end,3)];
LEOF2.R = [T1(end,7)  T1(end,8)  T1(end,9)];
MEO.R   = [T1(end,13) T1(end,14) T1(end,15)];
GEO.R   = [T1(end,19) T1(end,20) T1(end,21)];
SC.R    = [T1(end,25) T1(end,26) T1(end,27)];

LEOF1.V = [T1(end,4)  T1(end,5)  T1(end,6)];
LEOF2.V = [T1(end,10) T1(end,11) T1(end,12)];
MEO.V   = [T1(end,16) T1(end,17) T1(end,18)]; 
GEO.V   = [T1(end,22) T1(end,23) T1(end,24)];
SC.V    = [T1(end,28) T1(end,29) T1(end,30)];

% Circularize at perigee
DELTAV = DELTAV + abs(CIRCP(SC.rp) - CHELPP(SC.ra, SC.rp));

% RAAN change
ALPHA1 = ALPHAR(SC.i, abs(LEOF1.RAAN - LEOF2.RAAN));
DELTAV = DELTAV + DVRAAN(ALPHA1, norm(LEOF1.V));

% Update SC
SC.e    = 0;
SC.ra   = SC.rp;
SC.h    = sqrt(MU*SC.rp);
SC.RAAN = LEOF2.RAAN;

% Temp Circular Orbit that has radius equivalent to ra of LEO 2
TORB.rp = LEOF2.ra;
TORB.ra = LEOF2.ra;
TORB.h  = sqrt(MU*LEOF2.ra);

% Hohmann
[dv,dt,ecc,h,a,T,rp,ra] = Hohmann(SC, TORB);
DELTAV  = DELTAV + dv;

% Spacecraft Parameters
CTSCa   = COES(SC.R,SC.V);
SC.T    = T;          % s
SC.ra   = ra;         % km
SC.rp   = rp;         % km
SC.a    = a;          % km
SC.e    = ecc;               
SC.h    = h;          % km^2/s
SC.w    = LEOF1.w;    % rads
SC.Me   = CTSCa(12);  % rads

% Initial Eccentric Anomaly Guess
E0SC  = EGuess(SC.Me, SC.e);

% Final Eccentric Anomaly
ESC  = newtonsMethod(F, FPRIME, E0SC, TOL, SC.Me, SC.e);

% Get True Anomaly
SC.TA = THETA(ESC, SC.e); % rads

% Perifocal Vectors
R0SC   = RX(SC.h, SC.e, SC.TA); % km
V0SC   = VX(SC.h, SC.e, SC.TA); % km/s

% ECI State Vectors
Qsc    = QXX(SC.w, SC.i, SC.RAAN);
SC.R   = R_ECI(Qsc, R0SC); % km
SC.V   = v_ECI(Qsc, V0SC); % km/s

% Hohman-Transfer prop
timespan = [0 dt]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R SC.V];

% Periapse/apoapse prop
[tT1a, T1a] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
% **SC now on circularized LEO 2 with ra radius 
LEOF1.R = [T1a(end,1)  T1a(end,2)  T1a(end,3)];
LEOF2.R = [T1a(end,7)  T1a(end,8)  T1a(end,9)];
MEO.R   = [T1a(end,13) T1a(end,14) T1a(end,15)];
GEO.R   = [T1a(end,19) T1a(end,20) T1a(end,21)];
SC.R    = [T1a(end,25) T1a(end,26) T1a(end,27)];

LEOF1.V = [T1a(end,4)  T1a(end,5)  T1a(end,6)];
LEOF2.V = [T1a(end,10) T1a(end,11) T1a(end,12)];
MEO.V   = [T1a(end,16) T1a(end,17) T1a(end,18)]; 
GEO.V   = [T1a(end,22) T1a(end,23) T1a(end,24)];
SC.V    = [T1a(end,28) T1a(end,29) T1a(end,30)];

% Prop on circularized orbit back to LEO 2 apogee, tranfer to LEO 2 orbit
CTSCT    = COES(SC.R, SC.V);
CTL2T    = COES(LEOF2.R, LEOF2.V);
timespan = [0 (CTSCT(7) - (CTL2T(8) - CTSCT(8)) )];

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% Periapse/apoapse prop
[tT1b, T1b] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
% **SC now on LEO 2
LEOF1.R = [T1b(end,1)  T1b(end,2)  T1b(end,3)];
LEOF2.R = [T1b(end,7)  T1b(end,8)  T1b(end,9)];
MEO.R   = [T1b(end,13) T1b(end,14) T1b(end,15)];
GEO.R   = [T1b(end,19) T1b(end,20) T1b(end,21)];
SC.R    = [T1b(end,25) T1b(end,26) T1b(end,27)];

LEOF1.V = [T1b(end,4)  T1b(end,5)  T1b(end,6)];
LEOF2.V = [T1b(end,10) T1b(end,11) T1b(end,12)];
MEO.V   = [T1b(end,16) T1b(end,17) T1b(end,18)]; 
GEO.V   = [T1b(end,22) T1b(end,23) T1b(end,24)];
SC.V    = [T1b(end,28) T1b(end,29) T1b(end,30)];

% delta V for transfer
DELTAV = DELTAV + abs(LEOF2.h/LEOF2.ra - SC.h/SC.ra);


%% LEO 1 - LEO 2 Transfer Animation

% clear SA D vid
% SA = ones(length(T1(:,1))+length(T1a(:,1))+length(T1b(:,1)),30);
% for i = 1:30
%     tmp = {[T1(:,i)],[T1a(:,i)],[T1b(:,i)]};
%     SA(:,i) = vertcat(tmp{:});
% end
% 
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve5 = animatedline('LineWidth',0.25,'Color','y');
% 
% set(gca,'XLim',[min([min(SA(:,1)) min(SA(:,7)) min(SA(:,25))])...
%                 max([max(SA(:,1)) max(SA(:,7)) max(SA(:,25))])],...
%         'YLim',[min([min(SA(:,2)) min(SA(:,8)) min(SA(:,26))])...
%                 max([max(SA(:,2)) max(SA(:,8)) max(SA(:,26))])],...
%         'ZLim',[min([min(SA(:,3)) min(SA(:,9)) min(SA(:,27))])...
%                 max([max(SA(:,3)) max(SA(:,9)) max(SA(:,27))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(432, -75);
% hold on
% 
% 
% for i = 1:5:length(SA(:,1))
%     addpoints(curve1,SA(i,1), SA(i,2), SA(i,3));
%     addpoints(curve2,SA(i,7), SA(i,8), SA(i,9));
%     addpoints(curve5,SA(i,25),SA(i,26),SA(i,27));
%     
%     head1   = scatter3(SA(i,1),  SA(i,2),  SA(i,3),'MarkerEdgeColor','b');
%     head2   = scatter3(SA(i,7),  SA(i,8),  SA(i,9),'MarkerFaceColor', 'r','MarkerEdgeColor','r');
%     spccrft = scatter3(SA(i,25), SA(i,26), SA(i,27),'w>');
%     legend('LEO_1','LEO_2','S/C_T_R_A_N_S_F_E_R',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'TRANSTAGE_4_ _D_E_B','S/C')
%     legend('TextColor','white')
%     title('LEO 1 - LEO 2 Transfer')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(spccrft);
% end
% 
% vid = VideoWriter('LEO1_to_LEO2.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:5:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


%% Phasing Part 1 - Set Up

% Update Spacecraft to LEO 2 parameters
CTSC    = COES(SC.R, SC.V);
CTL2    = COES(LEOF2.R, LEOF2.V);
SC.T    = LEOF2.T;    % s
SC.i    = LEOF2.i;    % rads
SC.ra   = LEOF2.ra;   % km
SC.rp   = LEOF2.rp;   % km
SC.a    = LEOF2.a;    % km
SC.e    = LEOF2.e;               
SC.h    = LEOF2.h;    % km^2/s
SC.RAAN = LEOF2.RAAN; % rads
SC.w    = LEOF2.w;    % rads
SC.Me   = CTSC(12);   % rads
SC.TA   = pi + (CTL2(6) - CTSC(6));

% Perifocal Vectors
R0SC  = RX(SC.h, SC.e, SC.TA); % km
V0SC  = VX(SC.h, SC.e, SC.TA); % km/s

% ECI State Vectors
QSC   = QXX(SC.w, SC.i, SC.RAAN);
SC.R  = R_ECI(QSC, R0SC);     % km
SC.V  = v_ECI(QSC, V0SC);     % km/s

% Propagate to perigee
CTSCa = COES(SC.R, SC.V);
timespan = [0 (SC.T-CTSCa(8))]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R SC.V];

% Periapse/apoapse prop
[tPH1, PH1] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [PH1(end,1)  PH1(end,2)  PH1(end,3)];
LEOF2.R = [PH1(end,7)  PH1(end,8)  PH1(end,9)];
MEO.R   = [PH1(end,13) PH1(end,14) PH1(end,15)];
GEO.R   = [PH1(end,19) PH1(end,20) PH1(end,21)];
SC.R    = [PH1(end,25) PH1(end,26) PH1(end,27)];

LEOF1.V = [PH1(end,4)  PH1(end,5)  PH1(end,6)];
LEOF2.V = [PH1(end,10) PH1(end,11) PH1(end,12)];
MEO.V   = [PH1(end,16) PH1(end,17) PH1(end,18)]; 
GEO.V   = [PH1(end,22) PH1(end,23) PH1(end,24)];
SC.V    = [PH1(end,28) PH1(end,29) PH1(end,30)];


%% Phasing Part 2 - Catch LEO 2 Debris

% **INNER/OUTER TRANSFER
CTSCb = COES(SC.R, SC.V);
CTL2b = COES(LEOF2.R, LEOF2.V);

timespan = [0 (LEOF2.T - CTL2b(8) + LEOF2.T)];

% Update Spacecraft to Transfer Trajectory Parameters
SC.T    = timespan(2);       % s
SC.rp   = LEOF2.rp;          % km
SC.a    = A2(SC.T);          % km
SC.ra   = RA(SC.a, SC.rp);   % km
SC.e    = ECC(SC.ra, SC.rp);               
SC.h    = H(SC.e, SC.rp);    % km^2/s
SC.Me   = CTSCb(12);         % rads

% Delta-V for Phase
V1 = LEOF2.h/LEOF2.rp;
V2 = SC.h/SC.rp;
DELTAV = DELTAV + 2*abs(V1 - V2);

% Initial Eccentric Anomaly Guess
E0SC  = EGuess(SC.Me, SC.e);

% Final Eccentric Anomaly
ESC1  = newtonsMethod(F, FPRIME, E0SC, TOL, SC.Me, SC.e);

% Get True Anomaly
SC.TA = THETA(ESC1, SC.e);     % rads

% Perifocal Vectors
R0SC  = RX(SC.h, SC.e, SC.TA); % km
V0SC  = VX(SC.h, SC.e, SC.TA); % km/s

% ECI State Vectors
QSC   = QXX(SC.w, SC.i, SC.RAAN);
SC.R  = R_ECI(QSC, R0SC);     % km
SC.V  = v_ECI(QSC, V0SC);     % km/s

% Propagate Phase Maneuver
% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R SC.V];

% Periapse/apoapse prop
[tPH2, PH2] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [PH2(end,1)  PH2(end,2)  PH2(end,3)];
LEOF2.R = [PH2(end,7)  PH2(end,8)  PH2(end,9)];
MEO.R   = [PH2(end,13) PH2(end,14) PH2(end,15)];
GEO.R   = [PH2(end,19) PH2(end,20) PH2(end,21)];
SC.R    = [PH2(end,25) PH2(end,26) PH2(end,27)];

LEOF1.V = [PH2(end,4)  PH2(end,5)  PH2(end,6)];
LEOF2.V = [PH2(end,10) PH2(end,11) PH2(end,12)];
MEO.V   = [PH2(end,16) PH2(end,17) PH2(end,18)]; 
GEO.V   = [PH2(end,22) PH2(end,23) PH2(end,24)];
SC.V    = [PH2(end,28) PH2(end,29) PH2(end,30)];


%% Periapse Prop + Phasing Maneuver Animation

% clear SA D vid
% SA = ones(length(PH1(:,1))+length(PH2(:,1)),30);
% for i = 1:30
%     tmp = {[PH1(:,i)],[PH2(:,i)]};
%     SA(:,i) = vertcat(tmp{:});
% end
% 
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve5 = animatedline('LineWidth',0.25,'Color','y');
% 
% set(gca,'XLim',[min([min(SA(:,1)) min(SA(:,7)) min(SA(:,25))])...
%                 max([max(SA(:,1)) max(SA(:,7)) max(SA(:,25))])],...
%         'YLim',[min([min(SA(:,2)) min(SA(:,8)) min(SA(:,26))])...
%                 max([max(SA(:,2)) max(SA(:,8)) max(SA(:,26))])],...
%         'ZLim',[min([min(SA(:,3)) min(SA(:,9)) min(SA(:,27))])...
%                 max([max(SA(:,3)) max(SA(:,9)) max(SA(:,27))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(534, 68);
% hold on
% 
% 
% for i = 1:5:length(SA(:,1))
%     addpoints(curve1,SA(i,1), SA(i,2), SA(i,3));
%     addpoints(curve2,SA(i,7), SA(i,8), SA(i,9));
%     addpoints(curve5,SA(i,25),SA(i,26),SA(i,27));
%     
%     head1   = scatter3(SA(i,1),  SA(i,2),  SA(i,3),'MarkerEdgeColor','b');
%     head2   = scatter3(SA(i,7),  SA(i,8),  SA(i,9),'MarkerFaceColor', 'r','MarkerEdgeColor','r');
%     spccrft = scatter3(SA(i,25), SA(i,26), SA(i,27),'w>');
%     legend('LEO_1','LEO_2','S/C_P_H_A_S_I_N_G',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'TRANSTAGE_4_ _D_E_B','S/C')
%     legend('TextColor','white')
%     title('LEO 2 Phase')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(spccrft);
% end
% 
% vid = VideoWriter('LEO2_PHASE.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:5:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


%% 2nd LEO Object - 5 PERIODS

% Set LEO 2 & SC Equivalent
SC.R = LEOF2.R; 
SC.V = LEOF2.V;

timespan = [0 LEOF2.T*5]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% LEOF2 prop
[t2, S2] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [S2(end,1)  S2(end,2)  S2(end,3)];
LEOF2.R = [S2(end,7)  S2(end,8)  S2(end,9)];
MEO.R   = [S2(end,13) S2(end,14) S2(end,15)];
GEO.R   = [S2(end,19) S2(end,20) S2(end,21)];
SC.R    = [S2(end,25) S2(end,26) S2(end,27)];

LEOF1.V = [S2(end,4)  S2(end,5)  S2(end,6)];
LEOF2.V = [S2(end,10) S2(end,11) S2(end,12)];
MEO.V   = [S2(end,16) S2(end,17) S2(end,18)]; 
GEO.V   = [S2(end,22) S2(end,23) S2(end,24)];
SC.V    = [S2(end,28) S2(end,29) S2(end,30)];


%% LEO 2 Animation

% clear SA D vid
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% 
% set(gca,'XLim',[min([min(S2(:,1)) min(S2(:,7)) min(S2(:,13)) min(S2(:,19))])...
%                 max([max(S2(:,1)) max(S2(:,7)) max(S2(:,13)) max(S2(:,19))])],...
%         'YLim',[min([min(S2(:,2)) min(S2(:,8)) min(S2(:,14)) min(S2(:,20))])...
%                 max([max(S2(:,2)) max(S2(:,8)) max(S2(:,14)) max(S2(:,20))])],...
%         'ZLim',[min([min(S2(:,3)) min(S2(:,9)) min(S2(:,15)) min(S2(:,21))])...
%                 max([max(S2(:,3)) max(S2(:,9)) max(S2(:,15)) max(S2(:,21))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(307, 60);
% % view(0, 90);
% hold on
% 
% 
% for i = 1:6:length(S2(:,1))
%     addpoints(curve1,S2(i,1), S2(i,2), S2(i,3));
%     addpoints(curve2,S2(i,7), S2(i,8), S2(i,9));
%     addpoints(curve3,S2(i,13),S2(i,14),S2(i,15));
%     addpoints(curve4,S2(i,19),S2(i,20),S2(i,21));
%     
%     head1   = scatter3(S2(i,1),  S2(i,2),  S2(i,3), 'MarkerEdgeColor','b');
%     spccrft = scatter3(S2(i,7),  S2(i,8),  S2(i,9),'w>');
%     head3   = scatter3(S2(i,13), S2(i,14), S2(i,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c');
%     head4   = scatter3(S2(i,19), S2(i,20), S2(i,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m');
%     legend('LEO_1','LEO_2','MEO','GEO','DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'S/C & TRANSTAGE_4_ _D_E_B',...
%            'SATCOM_3_ _D_E_B','FENGYUN_2_H_ _D_E_B')
%     legend('TextColor','white')
%     title('LEO 2')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(spccrft);
%     delete(head3);
%     delete(head4);
% end
% 
% vid = VideoWriter('LEO2.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:6:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


% plot3(S2(:,1),  S2(:,2),  S2(:,3), 'b',...
%       S2(:,7),  S2(:,8),  S2(:,9), 'r',...
%       S2(:,13), S2(:,14), S2(:,15),'c',...
%       S2(:,19), S2(:,20), S2(:,21),'m')
% hold on
% scatter3(S2(end,1),  S2(end,2),  S2(end,3), 'MarkerFaceColor', 'b','MarkerEdgeColor','b')
% scatter3(S2(end,7),  S2(end,8),  S2(end,9), 'MarkerFaceColor', 'r','MarkerEdgeColor','r')
% scatter3(S2(end,13), S2(end,14), S2(end,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c')
% scatter3(S2(end,19), S2(end,20), S2(end,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m')
% scatter3(S2(end,25), S2(end,26), S2(end,27),'k>')


%% 2nd Transfer - Lamberts

% Prop MEO Object to its periapse
CTM1     = COES(MEO.R, MEO.V);
timespan = [0 (MEO.T-CTM1(8))]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% LEOF2 prop
[tT2, ST2] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [ST2(end,1)  ST2(end,2)  ST2(end,3)];
LEOF2.R = [ST2(end,7)  ST2(end,8)  ST2(end,9)];
MEO.R   = [ST2(end,13) ST2(end,14) ST2(end,15)];
GEO.R   = [ST2(end,19) ST2(end,20) ST2(end,21)];
SC.R    = [ST2(end,25) ST2(end,26) ST2(end,27)];

LEOF1.V = [ST2(end,4)  ST2(end,5)  ST2(end,6)];
LEOF2.V = [ST2(end,10) ST2(end,11) ST2(end,12)];
MEO.V   = [ST2(end,16) ST2(end,17) ST2(end,18)]; 
GEO.V   = [ST2(end,22) ST2(end,23) ST2(end,24)];
SC.V    = [ST2(end,28) ST2(end,29) ST2(end,30)];

STRTV    = norm(LEOF2.V);
timespan = [0 3*LEOF2.T/5]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% LEOF2 prop
[tT2a, ST2a] = ode45(@propagate, timespan, state, options, MU);

% Get positon where Lambert's transfer orbit will meet MEO Object
MEO_OBJ   = [ST2a(end,13) ST2a(end,14) ST2a(end,15)];

% **SC now on MEO
[V1, V2] = Lamberts(LEOF2.R, MEO_OBJ, timespan(2));

% Prop All Orbits
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' LEOF2.R' V1'];
[tT2b, ST2b] = ode45(@propagate, timespan, state, options, MU);

LEOF1.R  = [ST2b(end,1)  ST2b(end,2)  ST2b(end,3)];
LEOF2.R  = [ST2b(end,7)  ST2b(end,8)  ST2b(end,9)];
MEO.R    = [ST2b(end,13) ST2b(end,14) ST2b(end,15)];
GEO.R    = [ST2b(end,19) ST2b(end,20) ST2b(end,21)];
SC.R     = [ST2b(end,25) ST2b(end,26) ST2b(end,27)];

LEOF1.V  = [ST2b(end,4)  ST2b(end,5)  ST2b(end,6)];
LEOF2.V  = [ST2b(end,10) ST2b(end,11) ST2b(end,12)];
MEO.V    = [ST2b(end,16) ST2b(end,17) ST2b(end,18)]; 
GEO.V    = [ST2b(end,22) ST2b(end,23) ST2b(end,24)];
SC.V     = [ST2b(end,28) ST2b(end,29) ST2b(end,30)];

DELTAV   = DELTAV + abs(STRTV - norm(V1)) + abs(norm(MEO.V) - norm(V2));


%% LEO 2 - MEO Transfer Animation

% clear SA D vid
% SA = ones(length(ST2(:,1))+length(ST2b(:,1)),30);
% for i = 1:30
%     tmp = {[ST2(:,i)],[ST2b(:,i)]};
%     SA(:,i) = vertcat(tmp{:});
% end
% 
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% curve5 = animatedline('LineWidth',0.25,'Color','y');
% 
% set(gca,'XLim',[min([min(SA(:,1)) min(SA(:,7)) min(SA(:,13)) min(SA(:,19)) min(SA(end,25))])...
%                 max([max(SA(:,1)) max(SA(:,7)) max(SA(:,13)) max(SA(:,19)) max(SA(end,25))])],...
%         'YLim',[min([min(SA(:,2)) min(SA(:,8)) min(SA(:,14)) min(SA(:,20)) min(SA(end,26))])...
%                 max([max(SA(:,2)) max(SA(:,8)) max(SA(:,14)) max(SA(:,20)) max(SA(end,26))])],...
%         'ZLim',[min([min(SA(:,3)) min(SA(:,9)) min(SA(:,15)) min(SA(:,21)) min(SA(end,27))])...
%                 max([max(SA(:,3)) max(SA(:,9)) max(SA(:,15)) max(SA(:,21)) max(SA(end,27))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(-9, 68);
% hold on
% 
% for i = 1:3:length(SA(:,1))
%     addpoints(curve1,SA(i,1), SA(i,2), SA(i,3));
%     addpoints(curve2,SA(i,7), SA(i,8), SA(i,9));
%     addpoints(curve3,SA(i,13),SA(i,14),SA(i,15));
%     addpoints(curve4,SA(i,19),SA(i,20),SA(i,21));
%     addpoints(curve5,SA(i,25),SA(i,26),SA(i,27));
%     
%     head1   = scatter3(SA(i,1),  SA(i,2),  SA(i,3), 'MarkerEdgeColor','b');
%     head2   = scatter3(SA(i,7),  SA(i,8),  SA(i,9), 'MarkerEdgeColor','r');
%     head3   = scatter3(SA(i,13), SA(i,14), SA(i,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c');
%     head4   = scatter3(SA(i,19), SA(i,20), SA(i,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m');
%     spccrft = scatter3(SA(i,25), SA(i,26), SA(i,27),'w>');
%     legend('LEO_1','LEO_2','MEO','GEO','Transfer','DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'SATCOM_3_ _D_E_B','FENGYUN_2_H_ _D_E_B','S/C')
%     title('LEO 2 - MEO Transfer')
%     legend('TextColor','white')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(head3);
%     delete(head4);
%     delete(spccrft);
% end
% 
% vid = VideoWriter('LEO2_to_MEO.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:6:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


% plot3(ST2b(:,1),  ST2b(:,2),  ST2b(:,3), 'b',...
%       ST2b(:,7),  ST2b(:,8),  ST2b(:,9), 'r',...
%       ST2b(:,13), ST2b(:,14), ST2b(:,15),'c',...
%       ST2b(:,19), ST2b(:,20), ST2b(:,21),'m',...
%       ST2b(:,25), ST2b(:,26), ST2b(:,27),'g')
% hold on
% scatter3(ST2b(end,1),  ST2b(end,2),  ST2b(end,3), 'MarkerFaceColor', 'b','MarkerEdgeColor','b')
% scatter3(ST2b(end,7),  ST2b(end,8),  ST2b(end,9), 'MarkerFaceColor', 'r','MarkerEdgeColor','r')
% scatter3(ST2b(end,13), ST2b(end,14), ST2b(end,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c')
% scatter3(ST2b(end,19), ST2b(end,20), ST2b(end,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m')
% scatter3(ST2b(end,25), ST2b(end,26), ST2b(end,27),'k>')


%% MEO Object - 5 PERIODS

timespan = [0 MEO.T*5]; % sec
SC.R     = MEO.R;
SC.V     = MEO.V;

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% MEO prop
[T3, S3] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [S3(end,1)  S3(end,2)  S3(end,3)];
LEOF2.R = [S3(end,7)  S3(end,8)  S3(end,9)];
MEO.R   = [S3(end,13) S3(end,14) S3(end,15)];
GEO.R   = [S3(end,19) S3(end,20) S3(end,21)];
SC.R    = [S3(end,25) S3(end,26) S3(end,27)];

LEOF1.V = [S3(end,4)  S3(end,5)  S3(end,6)];
LEOF2.V = [S3(end,10) S3(end,11) S3(end,12)];
MEO.V   = [S3(end,16) S3(end,17) S3(end,18)];
GEO.V   = [S3(end,22) S3(end,23) S3(end,24)];
SC.V    = [S3(end,28) S3(end,29) S3(end,30)];


%% MEO Animation

% clear SA D vid
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% 
% set(gca,'XLim',[min([min(S3(:,1)) min(S3(:,7)) min(S3(:,13)) min(S3(:,19))])...
%                 max([max(S3(:,1)) max(S3(:,7)) max(S3(:,13)) max(S3(:,19))])],...
%         'YLim',[min([min(S3(:,2)) min(S3(:,8)) min(S3(:,14)) min(S3(:,20))])...
%                 max([max(S3(:,2)) max(S3(:,8)) max(S3(:,14)) max(S3(:,20))])],...
%         'ZLim',[min([min(S3(:,3)) min(S3(:,9)) min(S3(:,15)) min(S3(:,21))])...
%                 max([max(S3(:,3)) max(S3(:,9)) max(S3(:,15)) max(S3(:,21))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(307, 60);
% % view(0, 90);
% hold on
% 
% 
% for i = 1:18:length(S3(:,1))
%     addpoints(curve1,S3(i,1), S3(i,2), S3(i,3));
%     addpoints(curve2,S3(i,7), S3(i,8), S3(i,9));
%     addpoints(curve3,S3(i,13),S3(i,14),S3(i,15));
%     addpoints(curve4,S3(i,19),S3(i,20),S3(i,21));
%     
%     head1   = scatter3(S3(i,1),  S3(i,2),  S3(i,3), 'MarkerEdgeColor','b');
%     head2   = scatter3(S3(i,7),  S3(i,8),  S3(i,9), 'MarkerEdgeColor','r');
%     spccrft = scatter3(S3(i,13),  S3(i,14),  S3(i,15),'w>');
%     head4   = scatter3(S3(i,19), S3(i,20), S3(i,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m');
%     legend('LEO_1','LEO_2','MEO','GEO','DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'S/C & SATCOM_3_ _D_E_B','FENGYUN_2_H_ _D_E_B')
%     legend('TextColor','white')
%     title('MEO')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(spccrft);
%     delete(head4);
% end
% 
% vid = VideoWriter('MEO.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:18:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


% plot3(S3(:,1),  S3(:,2),  S3(:,3), 'b',...
%       S3(:,7),  S3(:,8),  S3(:,9), 'r',...
%       S3(:,13), S3(:,14), S3(:,15),'c',...
%       S3(:,19), S3(:,20), S3(:,21),'m')
% hold on
% scatter3(S3(end,1),  S3(end,2),  S3(end,3), 'MarkerEdgeColor','b')
% scatter3(S3(end,7),  S3(end,8),  S3(end,9), 'MarkerEdgeColor','r')
% scatter3(S3(end,13), S3(end,14), S3(end,15),'MarkerFaceColor', 'c','MarkerEdgeColor','c')
% scatter3(S3(end,19), S3(end,20), S3(end,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m')
% scatter3(S3(end,25), S3(end,26), S3(end,27),'k>')


%% MEO - GEO Transfer Lamberts

timespan = [0 7*LEOF2.T/6]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% MEO prop
[tT3, ST3] = ode45(@propagate, timespan, state, options, MU);

% Get positon where Lambert's transfer orbit will meet GEO Object
GEO_OBJ   = [ST3(end,19) ST3(end,20) ST3(end,21)];
STRTV     = norm(MEO.V);

% **SC now on GEO
[V1, V2] = Lamberts(MEO.R, GEO_OBJ, timespan(2));

% Prop All Orbits
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' MEO.R' V1'];
[tT3a, ST3a] = ode45(@propagate, timespan, state, options, MU);

LEOF1.R  = [ST3a(end,1)  ST3a(end,2)  ST3a(end,3)];
LEOF2.R  = [ST3a(end,7)  ST3a(end,8)  ST3a(end,9)];
MEO.R    = [ST3a(end,13) ST3a(end,14) ST3a(end,15)];
GEO.R    = [ST3a(end,19) ST3a(end,20) ST3a(end,21)];
SC.R     = [ST3a(end,25) ST3a(end,26) ST3a(end,27)];

LEOF1.V  = [ST3a(end,4)  ST3a(end,5)  ST3a(end,6)];
LEOF2.V  = [ST3a(end,10) ST3a(end,11) ST3a(end,12)];
MEO.V    = [ST3a(end,16) ST3a(end,17) ST3a(end,18)]; 
GEO.V    = [ST3a(end,22) ST3a(end,23) ST3a(end,24)];
SC.V     = [ST3a(end,28) ST3a(end,29) ST3a(end,30)];

DELTAV = DELTAV + abs(STRTV - norm(V1)) + abs(norm(GEO.V) - norm(V2));


%% MEO - GEO Tranfer Animation

% clear SA D vid
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% curve5 = animatedline('LineWidth',0.25,'Color','y');
% 
% set(gca,'XLim',[min([min(ST3a(:,1)) min(ST3a(:,7)) min(ST3a(:,13)) min(ST3a(:,19)) min(ST3a(end,25))])...
%                 max([max(ST3a(:,1)) max(ST3a(:,7)) max(ST3a(:,13)) max(ST3a(:,19)) max(ST3a(end,25))])],...
%         'YLim',[min([min(ST3a(:,2)) min(ST3a(:,8)) min(ST3a(:,14)) min(ST3a(:,20)) min(ST3a(end,26))])...
%                 max([max(ST3a(:,2)) max(ST3a(:,8)) max(ST3a(:,14)) max(ST3a(:,20)) max(ST3a(end,26))])],...
%         'ZLim',[min([min(ST3a(:,3)) min(ST3a(:,9)) min(ST3a(:,15)) min(ST3a(:,21)) min(ST3a(end,27))])...
%                 max([max(ST3a(:,3)) max(ST3a(:,9)) max(ST3a(:,15)) max(ST3a(:,21)) max(ST3a(end,27))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(307, 60);
% hold on
% 
% for i = 1:3:length(ST3a(:,1))
%     addpoints(curve1,ST3a(i,1), ST3a(i,2), ST3a(i,3));
%     addpoints(curve2,ST3a(i,7), ST3a(i,8), ST3a(i,9));
%     addpoints(curve3,ST3a(i,13),ST3a(i,14),ST3a(i,15));
%     addpoints(curve4,ST3a(i,19),ST3a(i,20),ST3a(i,21));
%     addpoints(curve5,ST3a(i,25),ST3a(i,26),ST3a(i,27));
%     
%     head1   = scatter3(ST3a(i,1),  ST3a(i,2),  ST3a(i,3), 'MarkerEdgeColor','b');
%     head2   = scatter3(ST3a(i,7),  ST3a(i,8),  ST3a(i,9), 'MarkerEdgeColor','r');
%     head3   = scatter3(ST3a(i,13), ST3a(i,14), ST3a(i,15),'MarkerEdgeColor','c');
%     head4   = scatter3(ST3a(i,19), ST3a(i,20), ST3a(i,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m');
%     spccrft = scatter3(ST3a(i,25), ST3a(i,26), ST3a(i,27),'w>');
%     legend('LEO_1','LEO_2','MEO','GEO','Transfer','DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_S_A_T_C_O_M_ _3','FENGYUN_2_H_ _D_E_B','S/C')
%     title('MEO - GEO Transfer')
%     legend('TextColor','white')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(head3);
%     delete(head4);
%     delete(spccrft);
% end
% 
% vid = VideoWriter('MEO_to_GEO.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:3:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);


% plot3(ST3a(:,1),  ST3a(:,2),  ST3a(:,3), 'b',...
%       ST3a(:,7),  ST3a(:,8),  ST3a(:,9), 'r',...
%       ST3a(:,13), ST3a(:,14), ST3a(:,15),'c',...
%       ST3a(:,19), ST3a(:,20), ST3a(:,21),'m')
% hold on
% scatter3(ST3a(end,1),  ST3a(end,2),  ST3a(end,3), 'MarkerEdgeColor','b')
% scatter3(ST3a(end,7),  ST3a(end,8),  ST3a(end,9), 'MarkerEdgeColor','r')
% scatter3(ST3a(end,13), ST3a(end,14), ST3a(end,15),'MarkerEdgeColor','c')
% scatter3(ST3a(end,19), ST3a(end,20), ST3a(end,21),'MarkerFaceColor', 'm','MarkerEdgeColor','m')
% scatter3(ST3a(end,25), ST3a(end,26), ST3a(end,27),'k>')


%% GEO Object - 5 PERIODS

timespan = [0 GEO.T*5]; % sec

% State Vector
state = [LEOF1.R' LEOF1.V' LEOF2.R' LEOF2.V' MEO.R' MEO.V' GEO.R' GEO.V' SC.R' SC.V'];

% last prop
[T4, S4] = ode45(@propagate, timespan, state, options, MU);

% Update Position and Velocity Vectors
LEOF1.R = [S4(end,1)  S4(end,2)  S4(end,3)];
LEOF2.R = [S4(end,7)  S4(end,8)  S4(end,9)];
MEO.R   = [S4(end,13) S4(end,14) S4(end,15)];
GEO.R   = [S4(end,19) S4(end,20) S4(end,21)];

LEOF1.V = [S4(end,4)  S1(end,5)  S4(end,6)];
LEOF2.V = [S4(end,10) S1(end,11) S4(end,12)];
MEO.V   = [S4(end,16) S1(end,17) S4(end,18)];
GEO.V   = [S4(end,22) S1(end,23) S4(end,24)];


%% GEO Animation

% clear SA D vid
% curve1 = animatedline('LineWidth',0.25,'Color','b');
% curve2 = animatedline('LineWidth',0.25,'Color','r');
% curve3 = animatedline('LineWidth',0.25,'Color','c');
% curve4 = animatedline('LineWidth',0.25,'Color','m');
% 
% set(gca,'XLim',[min([min(S4(:,1)) min(S4(:,7)) min(S4(:,13)) min(S4(:,19))])...
%                 max([max(S4(:,1)) max(S4(:,7)) max(S4(:,13)) max(S4(:,19))])],...
%         'YLim',[min([min(S4(:,2)) min(S4(:,8)) min(S4(:,14)) min(S4(:,20))])...
%                 max([max(S4(:,2)) max(S4(:,8)) max(S4(:,14)) max(S4(:,20))])],...
%         'ZLim',[min([min(S4(:,3)) min(S4(:,9)) min(S4(:,15)) min(S4(:,21))])...
%                 max([max(S4(:,3)) max(S4(:,9)) max(S4(:,15)) max(S4(:,21))])]);
% set(gca,'Color','k')
% set(gcf, 'Position', get(0, 'Screensize'));
% view(307, 60);
% hold on
% 
% 
% for i = 1:21:length(S4(:,1))
%     addpoints(curve1,S4(i,1), S4(i,2), S4(i,3));
%     addpoints(curve2,S4(i,7), S4(i,8), S4(i,9));
%     addpoints(curve3,S4(i,13),S4(i,14),S4(i,15));
%     addpoints(curve4,S4(i,19),S4(i,20),S4(i,21));
%     
%     head1   = scatter3(S4(i,1),  S4(i,2),  S4(i,3), 'MarkerEdgeColor','b');
%     head2   = scatter3(S4(i,7),  S4(i,8),  S4(i,9), 'MarkerEdgeColor','r');
%     head3   = scatter3(S4(i,13), S4(i,14), S4(i,15),'MarkerEdgeColor','c');
%     spccrft = scatter3(S4(i,19), S4(i,20), S4(i,21),'w>');
%     legend('LEO_1','LEO_2','MEO','GEO','DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_T_R_A_N_S_T_A_G_E_ _4',...
%            'DEORBIT_S_A_T_C_O_M_ _3','S/C & FENGYUN_2_H_ _D_E_B')
%     legend('TextColor','white')
%     title('GEO')
%     drawnow
%     D(i) = getframe(gcf);
% 
%     delete(head1);
%     delete(head2);
%     delete(head3);
%     delete(spccrft);
% end
% 
% vid = VideoWriter('GEO.avi');
% vid.FrameRate = 20;
% open(vid);
% for i = 1:21:length(D)
%     writeVideo(vid,D(i));
% end
% close(vid);

