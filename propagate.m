function dstate = propagate(time, state, muearth)

% ========================================================
%
% Used in junction with ODE45 to propagate all orbits 
%   with the same time step.
% 
%   time    - unused passed parameter
%   state   - state vector
%   muearth - Earth's Gravitational Parameter
%
% ========================================================

    %% LEO 1
    L1x = state(1);
    L1y = state(2);
    L1z = state(3);

    L1dx = state(4);
    L1dy = state(5);
    L1dz = state(6);

    L1r = norm([L1x L1y L1z]);

    L1ddx = -muearth*L1x/L1r^3;
    L1ddy = -muearth*L1y/L1r^3;
    L1ddz = -muearth*L1z/L1r^3;
    
    %% LEO 2
    L2x   = state(7);
    L2y   = state(8);
    L2z   = state(9);

    L2dx  = state(10);
    L2dy  = state(11);
    L2dz  = state(12);

    L2r   = norm([L2x L2y L2z]);

    L2ddx = -muearth*L2x/L2r^3;
    L2ddy = -muearth*L2y/L2r^3;
    L2ddz = -muearth*L2z/L2r^3;
    
    %% MEO 
    Mx   = state(13);
    My   = state(14);
    Mz   = state(15);

    Mdx  = state(16);
    Mdy  = state(17);
    Mdz  = state(18);

    Mr   = norm([Mx My Mz]);

    Mddx = -muearth*Mx/Mr^3;
    Mddy = -muearth*My/Mr^3;
    Mddz = -muearth*Mz/Mr^3;
    
    %% GEO
    Gx   = state(19);
    Gy   = state(20);
    Gz   = state(21);

    Gdx  = state(22);
    Gdy  = state(23);
    Gdz  = state(24);

    Gr   = norm([Gx Gy Gz]);

    Gddx = -muearth*Gx/Gr^3;
    Gddy = -muearth*Gy/Gr^3;
    Gddz = -muearth*Gz/Gr^3;
    
    %% SC
    SCx   = state(25);
    SCy   = state(26);
    SCz   = state(27);

    SCdx  = state(28);
    SCdy  = state(29);
    SCdz  = state(30);

    SCr   = norm([SCx SCy SCz]);

    SCddx = -muearth*SCx/SCr^3;
    SCddy = -muearth*SCy/SCr^3;
    SCddz = -muearth*SCz/SCr^3;

    dstate = [L1dx; L1dy; L1dz; L1ddx; L1ddy; L1ddz;
              L2dx; L2dy; L2dz; L2ddx; L2ddy; L2ddz;
              Mdx;  Mdy;  Mdz;  Mddx;  Mddy;  Mddz;
              Gdx;  Gdy;  Gdz;  Gddx;  Gddy;  Gddz;
              SCdx; SCdy; SCdz; SCddx; SCddy; SCddz;];
end

