function guess = EGuess(Me, e)

% ========================================
%
% Initial Eccentric Anomaly Guess
%
%   Me - Mean Anomaly
%   e  - eccentricity
%
% ========================================

if Me < pi
    guess = Me + e/2;
else
    guess = Me - e/2;
end
end

