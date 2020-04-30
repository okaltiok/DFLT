function [model, PARAMS, rmse, error] = main(model, PARAMS,experiment,trial)
% main - excecutes one iteration of the system that consists of: 1)
% initialization, 2) localization and tracking, 3) smoothing, 4) parameter
% estimation
%
% [model, PARAMS, rmse, error] = main(model, PARAMS,experiment,trial)
%
%    model      - a struct that contains the model parameters
%    PARAMS     - a struct that contains experiment parameters
%    experiment - experiment ID = {1,...,6}
%    trial      - trial ID = {0,...,3}
% 
% Returns:
%    model  - a struct that contains the new model parameter estimates
%    PARAMS - a struct that contains the updated experiment parameters
%    rmse   - (1x2) root mean square error of the tracking and smoothing
%    error  - a struct that contains the Euclidean distance error of the
%             tracking and smoothing filters for each time instant.
%
% Example : 
%    model = []; PARAMS = []; experiment = 3; trial = 2;
%    [model, PARAMS, rmse, error] = main(model, PARAMS,experiment,trial)
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



% help variables 
PLOT_ON = true;     % set PLOT_ON false to disable visualization 
plot_handle = []; 
plot_time = 1;
processing_time = zeros(4,1);



%% Initialization of the system
tic

if isempty(PARAMS)
    % Load experiment data
    experiment_name = {'open_ch1','open_ch4','open_ch16','apt_ch1','apt_ch4','apt_ch16'};        
    PARAMS = load(strcat('experiments\',experiment_name{experiment}),'EXPERIMENT');
         
    % calculate initial estimates for mu using the first 120 seconds of data.
    PARAMS.EXPERIMENT.mu_0 = calibrate(PARAMS, ...
            PARAMS.EXPERIMENT.data(PARAMS.EXPERIMENT.idx(2,1)+1:PARAMS.EXPERIMENT.idx(1,2)-1,:),120);
     
    % determine start and end of the trial
    idx = PARAMS.EXPERIMENT.idx(:,trial+1);
    PARAMS.EXPERIMENT.data = PARAMS.EXPERIMENT.data(idx(1):idx(2),:);
    
    % power spectral density of process noise,  q = 1e-2 approx 1.8 m/s^2 acceleration
    PARAMS.EXPERIMENT.q = 1e-2;
end

% Get experiment parameters
NODES_NUMBER = PARAMS.EXPERIMENT.nodes_number;
CHANNEL_NUMBER = PARAMS.EXPERIMENT.channel_number;
CHANNELS = PARAMS.EXPERIMENT.channels;
data = PARAMS.EXPERIMENT.data; 

% define link indexes
PARAMS.EXPERIMENT.nonDiagIdx = ~diag(true(NODES_NUMBER,1));
PARAMS.EXPERIMENT.linkIndex = reshape(1:NODES_NUMBER*(NODES_NUMBER-1),NODES_NUMBER-1,NODES_NUMBER);

% Initialize model parameters 
if isempty(model)
    model.phi = -5.*ones(NODES_NUMBER*(NODES_NUMBER-1),CHANNEL_NUMBER);
    model.lambda = 0.04.*ones(NODES_NUMBER*(NODES_NUMBER-1),CHANNEL_NUMBER);
    model.mu = PARAMS.EXPERIMENT.mu_0; 
    model.sigma2 = ones(NODES_NUMBER*(NODES_NUMBER-1),CHANNEL_NUMBER);
end

% Initialize radio tomographic imaging
PARAMS = initialize_rti(PARAMS,model);

% Initialize arrays for storing data
KK = size(data,1);
state.XX = zeros(4,KK); 
state.XX([1 3],:) = data(:,24:25)';
state.MM = zeros(4,KK);
state.PP = zeros(4,4,KK);
state.ZZ = zeros(2,KK/NODES_NUMBER);
state.RR = zeros(2,2,KK/NODES_NUMBER);
signal.RSS = zeros(NODES_NUMBER-1,KK);
signal.t = cumsum(data(:,1)');

% initialize state using the first RTI position estimate
measurement.ch = CHANNELS(data(NODES_NUMBER,3));               % frequency channel index
measurement.RSS = data(1:NODES_NUMBER,4:23)';                  % RSS measurement index
[~, measurement] = rti(measurement, model, PARAMS);

m = [measurement.z_rti(1) 0 measurement.z_rti(2) 0]';
P = kron(measurement.R_rti,[1 0; 0 0]) + kron(eye(2),[0 0; 0 1]);

% initialize measurement selection unit
filter_state = 2;
measurement.z_prev = measurement.z_rti; 
measurement.R_prev = measurement.R_rti;

processing_time(1) = toc;



%% run localization and tracking system

for k = NODES_NUMBER:NODES_NUMBER:KK
    tic
    
    % extract measurements of one communication cycle from data array
    measurement.dt = data(k-NODES_NUMBER+1:k,1);        % time index
    measurement.ch = CHANNELS(data(k,3));               % frequency channel index
    measurement.RSS = data(k-NODES_NUMBER+1:k,4:23)';   % RSS measurement index
    measurement.x = data(k-NODES_NUMBER+1:k,24:25);     % Position index

    % tracking filter
    [m,P,MM,PP] = tracking_filter(m,P,measurement,model,filter_state,PARAMS);

    % radio tomographic imaging
    [z_i, measurement] = rti(measurement, model, PARAMS);
        
    % measurement selection
    filter_state =  measurement_selection(m,P,measurement);
    
    % store measurements for post-processing
    signal.RSS(:,k-NODES_NUMBER+1:k) = reshape(measurement.RSS(~logical(eye(NODES_NUMBER))),NODES_NUMBER-1,NODES_NUMBER);
    state.MM(:,k-NODES_NUMBER+1:k) = MM;
    state.PP(:,:,k-NODES_NUMBER+1:k) = PP;
    state.ZZ(:,k/NODES_NUMBER) = measurement.z_rti;
    state.RR(:,:,k/NODES_NUMBER) = chol(measurement.R_rti);
    measurement.z_prev = measurement.z_rti; 
    measurement.R_prev = measurement.R_rti;
    
    processing_time(2) =  processing_time(2) + toc;
    
    % plot estimates every 0.5 seconds
    plot_time = plot_time + sum(measurement.dt);
    if PLOT_ON && plot_time > 0.5
        plot_handle = visualize(m,P,measurement,z_i,PARAMS,plot_handle);
        plot_time = 0;
    end
end



%% Smoothing and parameter estimation

tic
% Smoothen state estimates
state = rtss(state,data(:,1),PARAMS.EXPERIMENT.q); 
processing_time(3) = toc;

tic
% estimate mu, phi and sigma2 using the EM algorithm
model = em(model,state,signal,PARAMS);
processing_time(4) = toc;

% evaluate system performance
[rmse,error,PARAMS] = evaluate_performance(PLOT_ON,[experiment trial],state,signal,processing_time,PARAMS);

end




  