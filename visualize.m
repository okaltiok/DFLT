% visualize - illustrate the estimated image and location with respect to
% the true location (if known)
%
% handle = visualize(m,P,measurement,z_i,PARAMS,handle) 
%
%    m              - (4x1) state vector
%    P              - (4x4) state covariance
%    measurement    - a struct that contains the measurements
%    z_i            - (Nx1) RTI image
%    PARAMS         - a struct that contains experiment parameters
%    handle         - a handle to the chart line objects

% 
% Returns:
%    handle         - the updated handle to the chart line objects
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

function handle = visualize(m,P,measurement,z_i,PARAMS,handle) 
    x = measurement.x(end,:)';  % true position, x = nan(2,1) if it isn't known
    R = measurement.R_rti;      % covariance of RTI pos. est.
    z = measurement.z_rti;      % RTI pos. est.

    % form uncertainty ellipses
    theta = linspace(0,2*pi,100);
    rti_posterior = repmat(z,1,length(theta)) + 3.*chol(R,'lower')*[cos(theta);sin(theta)];
    ekf_posterior = repmat(m([1 3]),1,length(theta)) + 3.*chol(P([1 3],[1 3]),'lower')*[cos(theta);sin(theta)];

    try
        % update plot handles
        set(handle(1),'CData',reshape(z_i, length(PARAMS.RTI.yVals), length(PARAMS.RTI.xVals)));
        set(handle(2),'Xdata',ekf_posterior(1,:),'Ydata',ekf_posterior(2,:))
        set(handle(3),'Xdata',rti_posterior(1,:),'Ydata',rti_posterior(2,:))
        set(handle(4),'XData',m(1),'YData',m(3));
        set(handle(5),'XData',x(1),'YData',x(2));
    catch
        % get node positions
        nodeLocs = PARAMS.EXPERIMENT.nodeLocs;
        
        % set figure and axes
        fhandle = figure(1); 
        clf; axis xy; hold on; box on; 
        set(fhandle,'color',[1 1 1]);
        colormap jet;

        axis([min(nodeLocs(1,:))-0.5 max(nodeLocs(1,:))+0.5 min(nodeLocs(2,:))-0.5 max(nodeLocs(2,:))+0.5])
        set(gca,'TickLabelInterpreter','latex','fontsize',16);
        xlabel('x [m]','interpreter','latex');
        ylabel('y [m]','interpreter','latex');
        
        % plot RTI image on the background
        handle(1) = imagesc(PARAMS.RTI.xVals, PARAMS.RTI.yVals ,reshape(z_i, length(PARAMS.RTI.yVals), length(PARAMS.RTI.xVals))');
        
        % plot apartment layout on top of RTI image
        if strncmp(PARAMS.EXPERIMENT.layout,'apt_',4)
            layout = imread('experiments\apartment.jpg');
            layout(1:40,:,:) = [];
            [Ny,Nx,~] = size(layout);
            pixel_size =  0.0184;
            
            x = (1:Nx)*pixel_size;
            y = (1:Ny)*pixel_size;
            x = x - 18*pixel_size;
            y = y - 16*pixel_size;
            y = fliplr(y);
            
            im = imagesc(x,y,layout);
            im.AlphaData = 0.5;
        end
        
        % plot nodes and reference positions
        idx = find(diff(PARAMS.EXPERIMENT.data(:,24)) == 0 & diff(PARAMS.EXPERIMENT.data(:,25)) == 0);
        referencePos = unique(PARAMS.EXPERIMENT.data(idx,24:25),'rows')';
        plot(referencePos(1,:),referencePos(2,:),'k+','markersize',8,'linewidth',1.5);
        plot(nodeLocs(1,:),nodeLocs(2,:),'ks','markersize',8,'linewidth',2)

        % estimates and ground truth position
        handle(2) = plot(ekf_posterior(1,:),ekf_posterior(2,:),'w','linewidth',1);
        handle(3) = plot(rti_posterior(1,:),rti_posterior(2,:),'w--','linewidth',2);
        handle(4) = plot(m(1),m(3),'wx','MarkerSize',16,'linewidth',2);
        handle(5) = plot(x(1),x(2),'wo','MarkerSize',16,'linewidth',2);
    end
    drawnow
end