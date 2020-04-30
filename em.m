% em - estimates the model parameters \theta = [\mu \phi \sigma]' of the 
% exponential model
%
% model = em(model,state,signal,PARAMS)
%
%    model          - a struct that contains the model parameters
%    state          - a struct that contains the state estimates
%    signal         - a struct that contains the measurements
%    PARAMS         - a struct that contains experiment parameters
% 
% Returns:
%    model          - the estimated model parameters
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

function model = em(model,state,signal,PARAMS)

% The G matrix in the EM-alogrithm can be close to singular so turn the 
% Matlab warning off
warning('off', 'MATLAB:singularMatrix')
warning('off', 'MATLAB:nearlySingularMatrix')

NODES_NUMBER = PARAMS.EXPERIMENT.nodes_number;
CHANNEL_NUMBER = PARAMS.EXPERIMENT.channel_number;
nodeLocs = PARAMS.EXPERIMENT.nodeLocs;
linkIndex = PARAMS.EXPERIMENT.linkIndex;
links = NODES_NUMBER*(NODES_NUMBER-1);

KK = size(signal.RSS,2);

% loop through all of the links
for ch = 1:CHANNEL_NUMBER
    for tx = 1:NODES_NUMBER
        % determine the sample indexes based on the transmitter and channel
        % identifiers
        idx = (ch-1)*NODES_NUMBER+tx:NODES_NUMBER*CHANNEL_NUMBER:KK;
        
        % Get the RSS and determine measurements that have not been received
        rss = signal.RSS(:,idx)';
        packet_received = rss ~= -128;
        
        % Get the smoothed state and covariance estimates
        m = state.MMS(:,idx);  
        P = state.PPS(:,:,idx); 
        
        p_tx = nodeLocs(:,tx);
        p_rx = nodeLocs; p_rx(:,tx) = [];
        
        % estimate the parameters one receiver at a time
        for i = 1:NODES_NUMBER-1
            l = linkIndex(i,tx);
            
            % calculate excess path length
            d_LoS = sqrt((p_rx(1,i) - p_tx(1)).^2 + (p_rx(2,i) - p_tx(2)).^2);
            d_tx = sqrt((p_tx(1)-m(1,packet_received(:,i))).^2 + (p_tx(2)-m(3,packet_received(:,i))).^2);
            d_rx = sqrt((p_rx(1,i)-m(1,packet_received(:,i))).^2 + (p_rx(2,i)-m(3,packet_received(:,i))).^2);
            delta = (d_tx + d_rx - d_LoS)';
                        
            % get the number of received measurements and the RSS
            n = sum(packet_received(:,i));
            y = rss(packet_received(:,i),i);
            
            % compute non-linear part of the exponential model
            h = exp(-delta./model.lambda(l,ch));
            
            % compute help variables used by the EM-algorithm
            D =  y'*y;
            B = [sum(y); h'*y];
            G = [n sum(h); sum(h) h'*h];
        
            delta_dx = - (p_tx(1)-m(1,packet_received(:,i)))./d_tx - (p_rx(1,i)-m(1,packet_received(:,i)))./d_rx;
            delta_dy = - (p_tx(2)-m(3,packet_received(:,i)))./d_tx - (p_rx(2,i)-m(3,packet_received(:,i)))./d_rx;
            h_dx = -1./model.lambda(l,ch).*h.*delta_dx';
            h_dy = -1./model.lambda(l,ch).*h.*delta_dy';
            k = 1;
            for j = 1:length(idx)
                if packet_received(j,i)
                    Hx = [h_dx(k) 0 h_dy(k) 0];
                    G(2,2) = G(2,2) + Hx*P(:,:,j)*Hx';
                    k = k + 1;
                end
            end
            
            % compute the parameter estimates
            theta_hat = G\B;
            sigma2_hat = (D - B'/G*B)/n;
            
            model.mu(l,ch) = theta_hat(1);
            model.sigma2(l,ch) = sigma2_hat;
            
            if min(delta) < model.lambda(l,ch)
                % update \phi only if the person has located within the
                % sensitivity region of the link
                model.phi(l,ch) = theta_hat(2);
            end
        end
    end

    % Apply shrinkage to the sample variance to assure non-zero
    % values and to pull the most extreme ones towards the mean value
    alpha = 0.1;    % shrinkage coefficient
    model.sigma2(:,ch) = (1-alpha)*model.sigma2(:,ch) + alpha*sum(model.sigma2(:,ch))/links;
end
