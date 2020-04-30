function mu_0 = calibrate(PARAMS,data,t_end)
    % calibrate - For each link, computes the mean of the RSS 
    % using data from the time window t = [0,...,t_end].
    %
    % mu_0 = calibrate(PARAMS)
    %
    %    PARAMS     - a struct that contains experiment parameters
    %    t_end      - a scalar defining length of the calibration window in
    %                 seconds
    % 
    % Returns:
    %    mu_0       - (S(S-1) x C) The mean RSS of each link
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

    
    if (nargin < 3)
        t_end = 120;    % by default, a 120 second time window is used
    end
    
    % Get experiment parameters
    S = PARAMS.EXPERIMENT.nodes_number;
    C = PARAMS.EXPERIMENT.channel_number;
    CHANNELS = PARAMS.EXPERIMENT.channels;

    % initialize arrays
    mu_0 = -128*ones(S*(S-1),C);
    linkIndex = reshape(1:S*(S-1),S-1,S);

    % extract first t_end seconds of data
    T = cumsum(data(:,1));
    T = T - T(1);
    t_idx = find(T>t_end,1,'first');
    if t_idx < S*C
        return
    elseif isempty(t_idx)
        t_idx = size(T,1);
    end
    data = data(1:t_idx,:);


    % Calculate mean of each link
    for k = 1:S*C
        tx = data(k,2);
        rss = data(k:S*C:end,4:23);
        ch = CHANNELS(data(k,3));

        RX = setdiff(1:S,tx);
        mu = -60*ones(S,1);

        for rx = 1:S
            if tx ~= rx
                idx = find(rss(:,rx) ~=-128);
                if ~isempty(idx)
                    mu(rx) = mean(rss(idx,rx));
                end
            end
        end     
        mu_0(linkIndex(:,tx),ch) = mu(RX);
    end
end





