% Description : 
%    This Matlab script is for evaluating the system performance. 
%
% Example : 
%    Select the section of code and press run (ctrl+enter)

% Author   : Ossi Kaltiokallio
%            Aalto University, School of Electrical Engineering
%            Department of Communications and Networking
%            Maarintie 8, 02150 Espoo
%            ossi.kaltiokallio@aalto.fi
% Last Rev : 5/2/2020
% Tested   : Matlab version 9.7.0.1190202 (R2019b)


% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.


%% 
%
% Run the system for N iterations with the given experiment and trial
%

N = 11;         % Iteration number
experiment = 3; % experiment = {1,...,6}
trial = 3;      % trial = {1,2,3} 
                % for trial=27, there is no ground truth trajectory

% By default, Matlab's multithreaded and single-threaded algorithms can
% yield different results. Regular Matlab uses BLAS for linear algebra 
% operations, and it is multithreaded. On the other hand, Matlab's
% multithreading runs single-threaded for each worker. To ensure the same
% results, run regular Matlab single threaded.
maxNumCompThreads(1);

model = []; 
params = [];
for n = 1:N
    [model, params, rmse] = main(model,params,experiment,trial);
end
% clear('experiment','trial','model','params','rmse','n','N')


%% 
%
% Run the system for N iterations with every experiment and trial. The code 
% below uses the parallel computing toolbox for faster excecution time. 
% You can also replace the parfor loop with a regular for-loop if your 
% Matlab doesn't have the required parallel computing toolbox. 
%

maxNumCompThreads(1);

idx = [reshape(repmat(1:6,3,1),18,1) reshape(repmat(1:3,1,6),18,1)];
L = size(idx,1);
N = 11;

result = zeros(L,2,N);
error = cell(L,N);
model = cell(L,1);
params = cell(L,1);
for n = 1:N
    fprintf('\nIteration: %d\n',n-1);
    parfor l = 1:L
        [model{l}, params{l}, rmse, e] = ...
            main(model{l}, params{l},idx(l,1),idx(l,2));
        result(l,:,n) = rmse;
        error{l,n} = e.m;
    end
end
clear('model','params','idx','l','L','n','N','rmse','e')


%% 
%
% Plot the RMSE as a function of iteration number and print a table that 
% summarizes the results in each experiment and trial. This section of code 
% is meant to be run after the section above is executed.
%

figure(1); clf; box on; grid on; hold on
plot(0:10,reshape(result(:,1,:),18,11),...
        'k','linewidth',0.5,'color',[0.7 0.7 0.7])
plot(0:10,(mean(reshape(result(:,1,:),18,11))),'k','linewidth',2)
set(gca,'TickLabelInterpreter','latex','fontsize',16);
ylabel('RMSE [m]','interpreter','latex');
xlabel('Iteration number','interpreter','latex');
legend('RMSE','Average','interpreter','latex')

for i = 1:18
    [~,idx(i)] = min(result(i,1,:));
end
I = reshape(1:18,3,6);

str = {'Ex1 (open)','Ex2 (open)','Ex3 (open)' ...
       'Ex4 (apt)', 'Ex5 (apt)', 'Ex6 (apt)'};
fprintf('%s %11s %15s %15s %15s\n','Experiment','Trial1','Trial2','Trial3','Average')
for j = 1:6
    e = [];
    fprintf([str{j},'\t'])
    for i = 1:3
        rmse = sqrt(mean((error{I(i,j),idx(I(i,j))}).^2));
        sigma = std(error{I(i,j),idx(I(i,j))});
        e = [e error{I(i,j),idx(I(i,j))}];
        fprintf(['%.2f ',char(177),' %.2f\t'],rmse*100,sigma*100)
    end
    rmse = sqrt(mean(e.^2));
    sigma = std(e);
    fprintf(['%.2f ',char(177),' %.2f\n'],rmse*100,sigma*100)
end
clear('i','I','j','e','rmse','sigma','idx','str')


% The code above should yield the following result, that is, the RMS
% position error ? standard deviation of the estimation error in
% centimeters.

% Experiment      Trial1          Trial2          Trial3         Average
% Ex1 (open)	41.51 ? 29.02	33.52 ? 19.64	31.67 ? 17.42	35.80 ? 22.61
% Ex2 (open)	19.76 ? 9.61	17.40 ? 8.17	22.41 ? 13.73	19.93 ? 10.78
% Ex3 (open)	17.37 ? 9.54	18.28 ? 9.46	17.72 ? 9.02	17.79 ? 9.36
% Ex4 (apt)     54.18 ? 31.15	48.16 ? 24.71	50.92 ? 26.34	51.20 ? 27.62
% Ex5 (apt)     35.55 ? 17.51	39.63 ? 21.70	36.41 ? 19.34	37.26 ? 19.65
% Ex6 (apt)     34.94 ? 21.17	35.53 ? 18.70	36.83 ? 19.63	35.78 ? 19.90



