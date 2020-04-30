% rtss - implements the backward smoothing recursion of the Rauch-Tung-
%        Striebel smoother
%
% state = rtss(state,dt,q)
%
%    state          - a struct that cointains the state estimates
%    dt             - a vector of sampling intervals
%    q              - power spectral density of process noise
% 
% Returns:
%    state          - updated struct that cointains the smoothed state
%                     estimates
%
% Author   : Ossi Kaltiokallio
%            Aalto University, School of Electrical Engineering
%            Department of Communications and Networking
%            Maarintie 8, 02150 Espoo
%            ossi.kaltiokallio@aalto.fi
% Last Rev : 6/2/2020
% Tested   : Matlab version 9.7.0.1190202 (R2019b)
%
% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

function state = rtss(state,dt,q)
    KK = size(dt,1);
    
    % initialize state and covariance
    state.MMS = state.MM;
    state.PPS = state.PP;
    m = state.MM(:,end);
    P = state.PP(:,:,end);

    % the backward recursion start at time sample K-1 and proceeds backward
    % in time to the first time instant
    for k =KK-1:-1:1
        
        % form transition and process noise matrices
        F = [1 dt(k) 0 0; ...
            0 1 0 0; ...
            0 0 1 dt(k); ...
            0 0 0 1];
        
        Q = q.*[dt(k)^3/3 dt(k)^2/2 0 0; ...
            dt(k)^2/2 dt(k) 0 0; ...
            0 0 dt(k)^3/3 dt(k)^2/2; ...
            0 0 dt(k)^2/2 dt(k)];

        % the RTSS algorithm
        m_ = F*state.MM(:,k);
        P_ = F*state.PP(:,:,k)*F' + Q;
        G = state.PP(:,:,k)*F'/P_;
        m = state.MM(:,k) + G*(m - m_);
        P = state.PP(:,:,k) + G*(P - P_)*G';

        state.MMS(:,k) = m;
        state.PPS(:,:,k) = P;
    end
end 