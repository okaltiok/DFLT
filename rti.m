% rti - estimate radio tomographic image from the RSS and
% estimate the person's location and covariance from the discretized image
%
% [image, measurement] = rti(measurement, model, PARAMS) 
%
%    measurement    - a struct that contains the measurements
%    model          - a struct that contains the model parameters
%    PARAMS         - a struct that contains experiment and imaging parameters
% 
% Returns:
%    image          - estimated RTI image
%    measurement    - measurement struct that is appended with position and
%                     covariance 
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

function [image, measurement] = rti(measurement, model, PARAMS) 
    gamma   = 0.70;                              % image threshold
    crlb = PARAMS.RTI.delta_p^2/12;               % CRLB of RTI
    
    % get RSS 
    Z = measurement.RSS(PARAMS.EXPERIMENT.nonDiagIdx);
    
    % Links that have received a packet
    idx = Z ~= -128;

    % create measurement vector and compute image
    z_r = sign(model.phi(idx,measurement.ch)) .* (Z(idx) - model.mu(idx,measurement.ch));
    image = PARAMS.RTI.projection(:,idx,measurement.ch) * z_r;
    
    % find pixels with intensity higher than threshold
    idx = image >= gamma*max(image);
    
    % calculate weight of pixels
    w = image(idx)./sum(image(idx));
    Y = PARAMS.RTI.pixelCoords(:,idx);
    n = size(Y,1);
    
    % weighted sample mean
    z_rti = zeros(n,1);
    for j = 1:n
        z_rti(j) = Y(j,:)*w;
    end
    
    % weighted sample covariance
    R_rti = zeros(n);
    for j = 1:n
        for k = 1:n
            R_rti(j,k) = w'.*(Y(j,:) - z_rti(j))*(Y(k,:) - z_rti(k))';
        end
    end
    
    % check that position estimate is a number
    if isnan(z_rti)
        [image, z_rti, R_rti] = default_est(PARAMS.RTI.pixels,PARAMS.RTI.pixelCoords);
    end
    
    % check that x- and y-coordinates of R_rti are larger than CRLB
    flag = false;
    [U,R,V] = svd(R_rti);
    if R(1,1) < crlb
        R(1,1) = crlb;
        flag = true;
    end
    if R(2,2) < crlb
        R(2,2) = crlb;
        flag = true;
    end
    if flag
        R_rti = U*R*V';
    end
    
    % check that R_rti is positive definite.
    [~,flag] = chol(R_rti,'lower');
    if flag ~= 0
        % Covariance matrix is not positive definite, use the CRLB value
        % instead. The Cramer-Rao lower bound (CRLB) of RTI is defined by
        % the pixel size (delta_p = 0.25) and it is delta_p^2/12 = 0.0052;
        R_rti = [crlb 0; 0 crlb];
    end
    
	% Append the position and covariance estimates to the measurement struct 
    measurement.z_rti = z_rti;
    measurement.R_rti = R_rti;
end


function [image, z_rti, R_rti] = default_est(n,pixelCoords)
    image = zeros(n,1);
    w = 1/n.*ones(n,1);
    z_rti = pixelCoords*w; 
    R_rti = 1e3*eye(2); 
end