% measurement_selection - determine the filter state based on the current
% state estimate and current and past RTI position estimates
%
% filter_state = measurement_selection(m,P,measurement) 
%
%    m              - (4x1) state vector
%    P              - (4x4) state covariance
%    measurement    - a struct that contains the measurements
% 
% Returns:
%    filter_state   - state of the filter {1,2,3}
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

function filter_state = measurement_selection(m,P,measurement)


threshold = 5.9915;         % threshold = chi2inv(.95,2)
H = [1 0 0 0; ...           % measurement model
     0 0 1 0];

% compute normalized innovation square (NIS), that is, square of the 
% Mahalanobis distance
epsilon_1 = (measurement.z_rti - H*m)'/ ...
            (H * P * H' + measurement.R_rti)*...
            (measurement.z_rti - H*m);

epsilon_2 = (measurement.z_rti - measurement.z_prev)'/...
            (measurement.R_prev + measurement.R_rti)*...
            (measurement.z_rti - measurement.z_prev);

        
% logic for determining filter state.
if epsilon_1 > threshold && epsilon_2 <= threshold
    
    % In this state, it is likely that the filter has diverged. In the next
    % filter recursion, only use RTI position estimates when updating the
    % filter.
    filter_state = 1;
    
elseif epsilon_1 <= threshold
    
    % Normal operation, use both RTI and RSS when updating the filter.
    filter_state = 2;
    
else
    
    % In this state, it is likely that the RTI estimate is inaccurate. In 
    % the next filter recursion, only use the RSS when updating the filter.
    filter_state = 3;
end