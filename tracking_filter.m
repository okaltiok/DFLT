% tracking_filter - perform prediction and update steps of the tracking
% filter for one communication cycle
%
% [m,P,MM,PP] = tracking_filter(m,P,measurement,model,filter_state,PARAMS)
%
%    m              - (4x1) state vector
%    P              - (4x4) state covariance
%    measurement    - a struct that contains the measurements
%    model          - a struct that contains the model parameters
%    filter_state   - a variable which indicates the model structure
%    PARAMS         - a struct that contains experiment parameters
% 
% Returns:
%    m              - (4x1) updated state vector.
%    P              - (4x4) updated state covariance
%    MM             - (4xN) an array for storing m at tx = 1,...,N
%    PP             - (4x4xN) an array for storing P at tx = 1,...,N
%
% Author   : Ossi Kaltiokallio
%            Aalto University, School of Electrical Engineering
%            Department of Communications and Networking
%            Maarintie 8, 02150 Espoo
%            ossi.kaltiokallio@aalto.fi
% Last Rev : 5/2/2020
% Tested   : Matlab version 9.7.0.1190202 (R2019b)
%
% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.

function [m,P,MM,PP] = tracking_filter(m,P,measurement,model,filter_state,PARAMS)
    q = PARAMS.EXPERIMENT.q;
    NODES_NUMBER = PARAMS.EXPERIMENT.nodes_number;
    nodeLocs = PARAMS.EXPERIMENT.nodeLocs;
    linkIndex = PARAMS.EXPERIMENT.linkIndex;

    % Read measurements
    z_rti = measurement.z_rti;
    R_rti = measurement.R_rti;
    dt = measurement.dt;
    ch = measurement.ch;
    RSS = measurement.RSS;
    
    % Initialize arrays to store state and covariance estimates
    MM = zeros(size(m,1),NODES_NUMBER);
    PP = zeros(size(m,1),size(m,1),NODES_NUMBER);
    
    % Initialize complete residual, measurement and covariance matrices
    HH = vertcat([1 0 0 0; 0 0 1 0], zeros(NODES_NUMBER-1,size(m,1)));
    RR = zeros(NODES_NUMBER+1,NODES_NUMBER+1);
    NU = zeros(NODES_NUMBER+1,1);
    RR_idx = boolean(blkdiag(zeros(2),eye(NODES_NUMBER-1)));
    meas_idx = true(NODES_NUMBER+1,1);
    
    for tx = 1:NODES_NUMBER
        rss =  RSS(:,tx);
        rss(tx) = [];
        packet_received = rss ~= -128;
           
        if dt(tx) <= 0, dt(tx) = 1e-9; end
        
        % form transition and process noise matrices
        F = [1 dt(tx) 0 0; ...
            0 1 0 0; ...
            0 0 1 dt(tx); ...
            0 0 0 1];
        
        Q = q.*[dt(tx)^3/3 dt(tx)^2/2 0 0; ...
            dt(tx)^2/2 dt(tx) 0 0; ...
            0 0 dt(tx)^3/3 dt(tx)^2/2; ...
            0 0 dt(tx)^2/2 dt(tx)];
        
        % Prediction step of the filter
        m_ = F*m;
        P_ = F*P*F' + Q;
       
        % Update state estimates if measurements are received
        if sum(double(packet_received)) > 0
            if filter_state == 1
                % Use RTI pos. est. in filter update --> KF
                R = R_rti;
                H = [1 0 0 0; 0 0 1 0];
                nu = z_rti - H*m_;
            else
                % Use RSS and/or RTI pos. est. in filter update --> EKF
                
                % get TX and RXs coordinates
                p_tx = nodeLocs(:,tx);
                p_rx = nodeLocs; p_rx(:,tx) = [];
                
                % calculate excess path length
                d_tx = sqrt((m_(1)-p_tx(1)).^2 + (m_(3)-p_tx(2)).^2);
                d_rx = sqrt((m_(1)-p_rx(1,:)).^2 + (m_(3)-p_rx(2,:)).^2);
                d_LoS = sqrt((p_rx(1,:) - p_tx(1)).^2 + (p_rx(2,:) - p_tx(2)).^2);
                delta = (d_tx + d_rx - d_LoS)';
                
                % Form complete residual, measurement and covariance matrices
                h = model.phi(linkIndex(:,tx),ch).*exp(-delta./model.lambda(linkIndex(:,tx),ch));
                NU(1:2) = z_rti - m_([1 3]); NU(3:end) = rss - h - model.mu(linkIndex(:,tx),ch);
                HH(3:end,[1 3]) = (-h./model.lambda(linkIndex(:,tx),ch).*((m_([1 3]) - p_tx)./d_tx + (m_([1 3]) - p_rx)./d_rx)');
                RR(1:2,1:2) = R_rti; RR(RR_idx) = model.sigma2(linkIndex(:,tx),ch);
                
                % extract sub model structure according to the filter state
                % and measured RSS
                if filter_state == 2
                    meas_idx(1:2) = true; 
                    meas_idx(3:end) = packet_received;
                elseif filter_state == 3
                    meas_idx(1:2) = false; 
                    meas_idx(3:end) = packet_received;
                end
                
                H = HH(meas_idx,:);
                R = RR(meas_idx,meas_idx);
                nu = NU(meas_idx,1);
            end
            
            % Update step of the filter
            S = H * P_ * H' + R;
            K = P_ *H' / S;
            m = m_ + K * nu;
            P = P_ - K * S * K';
  
            % check that P is positive definite.
            [~,flag] = chol(P,'lower');
            if flag ~= 0
                % try to force P to be positive definite and check again
                P = 0.5 .* (P + P');
                [~,flag] = chol(P,'lower');
                if flag ~= 0
                    m = m_; P = P_;
                end
            end
        else
            m = m_; 
            P = P_;
        end
        MM(:,tx) = m;
        PP(:,:,tx) = P;
    end
end
