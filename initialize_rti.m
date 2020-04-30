% initialize_rti - set up pixel locations and compute projection matrix
%
% PARAMS = initialize_rti(PARAMS,model)
%
%    model      - a struct that contains the model parameters
%    PARAMS     - a struct that contains experiment and imaging parameters
% 
% Returns:
%    PARAMS - a struct that contains the updated imaging parameters
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

function PARAMS = initialize_rti(PARAMS,model)
    NODES_NUMBER = PARAMS.EXPERIMENT.nodes_number;
    CHANNEL_NUMBER = PARAMS.EXPERIMENT.channel_number;
    nodeLocs = PARAMS.EXPERIMENT.nodeLocs;
    
    try
        % the variables have already been calculated 
        DistPixelAndNode = PARAMS.RTI.DistPixelAndNode;
        DistNodes = PARAMS.RTI.DistNodes;
        CovPixels = PARAMS.RTI.CovPixels;
        pixels = PARAMS.RTI.pixels;
    catch
        % image reconstruction parameters
        sigmax2       = 0.0005;       % variance of image value
        delta         = 0.5;          % pixel correlation distance
        delta_p       = 0.25;         % pixel size

        % Set up pixel locations
        xVals       = min(nodeLocs(1,:)) : delta_p : max(nodeLocs(1,:)) + (delta_p - mod(max(nodeLocs(1,:)),delta_p/2));
        yVals       = min(nodeLocs(2,:)) : delta_p : max(nodeLocs(2,:)) + (delta_p - mod(max(nodeLocs(2,:)),delta_p/2));
        [x, y] = meshgrid(xVals,yVals);
        pixelCoords = [x(:), y(:)]';
        pixels   = size(pixelCoords,2);
        
        % Find distances between pixels and transceivers
        DistPixels = L2_distance( pixelCoords, pixelCoords,1);
        DistPixelAndNode = L2_distance( pixelCoords, nodeLocs, 0);
        DistNodes = L2_distance( nodeLocs, nodeLocs, 1);
        CovPixels = (1/sigmax2).*inv(exp(-DistPixels./delta));
        
        PARAMS.RTI.sigmax2 = sigmax2;       
        PARAMS.RTI.delta = delta;         
        PARAMS.RTI.delta_p = delta_p;       
        
        PARAMS.RTI.pixelCoords = pixelCoords;
        PARAMS.RTI.pixels = pixels;
        PARAMS.RTI.xVals = xVals;
        PARAMS.RTI.yVals = yVals;
        PARAMS.RTI.DistNodes = DistNodes;
        PARAMS.RTI.DistPixelAndNode = DistPixelAndNode;
        PARAMS.RTI.CovPixels = CovPixels;
    end

    % Create mapping from link pair to link index, and backwards
    tx              = ones(NODES_NUMBER, 1)*(1:NODES_NUMBER);
    rx              = tx';
    idx             = ~diag(true(NODES_NUMBER,1));
    links           = NODES_NUMBER*(NODES_NUMBER-1);
    txForLinkNumber = tx(idx);
    rxForLinkNumber = rx(idx);
    
    % Compute projection matrix for each channel
    PARAMS.RTI.projection = zeros(pixels,links,CHANNEL_NUMBER);
    for ch = 1:CHANNEL_NUMBER
        W = zeros(links, pixels);
        % Calculate weight for each pixel and link
        for ln = 1:links
            txNum = txForLinkNumber(ln);
            rxNum = rxForLinkNumber(ln);
            d = DistPixelAndNode(:,txNum) + DistPixelAndNode(:,rxNum) - DistNodes(txNum,rxNum);
            W(ln, :) = exp(-d'./model.lambda(ln,ch));
        end
        
        % Measurement noise covariance
        R = diag(model.sigma2(:,ch));
        
        % Minimum mean square error estimate is
        PARAMS.RTI.projection(:,:,ch) = ((W'/R*W + CovPixels)\W')/R;
    end
end

