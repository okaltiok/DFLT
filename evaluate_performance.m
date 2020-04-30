% evaluate_performance - compute the RMS error of the filtered and smoothed
% state estimates
%
% [rmse,error] = evaluate_performance(PLOT_ON,ex,state,signal,processing_time,PARAMS) 
%
%    PLOT_ON            - if true, also plot the results
%    ex                 - experiment and trial ID
%    state              - a struct that cointains the state estimates
%    signal             - a struct that cointains the time vector
%    processing_time    - (3x1) vector of CPT times
%    PARAMS             - a struct that contains experiment and imaging parameters
% 
% Returns:
%    rmse         - RMS error of the filtered and smoothed state estimates
%    error        - RMS error of each time instant
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

function [rmse,error,PARAMS] = evaluate_performance(PLOT_ON,ex,state,signal,processing_time,PARAMS)
    
    % find elements with known ground truth position and calculate the RMSE
    pos_idx = find(~isnan(state.XX(1,:)));
    
    if ~isempty(pos_idx)
        error.m = sqrt(sum((state.MM([1 3],pos_idx) - state.XX([1 3],pos_idx)).^2,1));
        error.ms = sqrt(sum((state.MMS([1 3],pos_idx) - state.XX([1 3],pos_idx)).^2,1));
        rmse(1:2) = [sqrt(mean(error.m.^2)) sqrt(mean(error.ms.^2))];
    else
        rmse = nan(1,2);
        error = [];
    end
    
    fprintf('EX%dTR%d, RMSE filtered/smoothed: %.2f/%.2f [cm], CPU: %.2f %.2f %.2f %.2f [s]\n', ...
        ex(1),ex(2),rmse(1)*100,rmse(2)*100,processing_time)

    % plot estimates
    if PLOT_ON
        % Calculate positioning error of RTI, EKF and RTSS
        e1 = sqrt(sum((state.ZZ - state.XX([1 3],20:20:end)).^2,1));
        e2 = sqrt(sum((state.MM([1 3],:) - state.XX([1 3],:)).^2,1));
        e3 = sqrt(sum((state.MMS([1 3],:) - state.XX([1 3],:)).^2,1));
        
        signal.t(1) = 0;
        t = signal.t(1:20:end);
        
        try
            handle_x = PARAMS.VISUALIZE.handle_x;
            sigma2 = squeeze(state.RR(1,1,:))';

            set(handle_x(1),'Xlim',[signal.t(1) signal.t(end)],'Ylim',[floor(min(state.MM(1,:))) ceil(max(state.MM(1,:)))+1])
            set(handle_x(2), 'Ydata', [state.ZZ(1,:) - 3*sigma2 fliplr(state.ZZ(1,:) + 3*sigma2)])
            set(handle_x(4),'YData',state.MM(1,:)');
            set(handle_x(5),'YData',state.MMS(1,:)');
            
            handle_y = PARAMS.VISUALIZE.handle_y;
            sigma2 = squeeze(state.RR(2,2,:))';
             
            set(handle_y(1),'Xlim',[signal.t(1) signal.t(end)],'Ylim',[floor(min(state.MM(3,:))) ceil(max(state.MM(3,:)))+1])
            set(handle_y(2), 'Ydata', [state.ZZ(2,:) - 3*sigma2 fliplr(state.ZZ(2,:) + 3*sigma2)])
            set(handle_y(4),'YData',state.MM(3,:)');
            set(handle_y(5),'YData',state.MMS(3,:)');
            
            handle_e = PARAMS.VISUALIZE.handle_e;
            
            set(handle_e(1),'Xlim',[signal.t(1) signal.t(end)],'YScale','log','YLim',[0.01 10])
            set(handle_e(2),'YData',e1);
            set(handle_e(3),'YData',e2);
            set(handle_e(4),'YData',e3);
        catch
            figure(2); clf
            
            % x-coordinate
            handle_x(1) = subplot(311);
            sigma2 = squeeze(state.RR(1,1,:))';
            handle_x(2) = fill([t fliplr(t)],[state.ZZ(1,:) - 3*sigma2 fliplr(state.ZZ(1,:) + 3*sigma2)],'r','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
            hold on; box on; grid on
            handle_x(3) = plot(signal.t,state.XX(1,:)','k');
            handle_x(4) = plot(signal.t,state.MM(1,:)','b');
            handle_x(5) = plot(signal.t,state.MMS(1,:)','r:');
            set(gca,'Xlim',[signal.t(1) signal.t(end)],'Ylim',[floor(min(state.MM(1,:))) ceil(max(state.MM(1,:)))+1])
            set(gca,'TickLabelInterpreter','latex','fontsize',12);
            set(gca,'XtickLabel','')
            ylabel('x [m]','interpreter','latex');
            p = get(gca,'position'); set(gca,'position', [p(1:3) 0.27])
            legend('RTI $3\sigma$ CI','$x$','EKF','RTSS','interpreter','latex','orientation','horizontal')
            PARAMS.VISUALIZE.handle_x = handle_x;
            
            % y-coordinate
            handle_y(1) = subplot(312);
            sigma2 = squeeze(state.RR(2,2,:))';
            handle_y(2) = fill([t fliplr(t)],[state.ZZ(2,:) - 3*sigma2 fliplr(state.ZZ(2,:) + 3*sigma2)],'r','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
            hold on; box on; grid on
            handle_y(3) = plot(signal.t,state.XX(3,:)','k');
            handle_y(4) = plot(signal.t,state.MM(3,:)','b');
            handle_y(5) = plot(signal.t,state.MMS(3,:)','r:');
            set(gca,'Xlim',[signal.t(1) signal.t(end)],'Ylim',[floor(min(state.MM(3,:))) ceil(max(state.MM(3,:)))+1])
            set(gca,'TickLabelInterpreter','latex','fontsize',12);
            set(gca,'XtickLabel','')
            ylabel('y [m]','interpreter','latex');
            p = get(gca,'position'); set(gca,'position', [p(1:3) 0.27])
            legend('RTI $3\sigma$ CI','$y$','EKF','RTSS','interpreter','latex','orientation','horizontal')
            PARAMS.VISUALIZE.handle_y = handle_y;
            
            % positioning error
            handle_e(1) = subplot(313);
            handle_e(2) = plot(t,e1,'color',[1 0 0 0.25]);
            hold on; box on; grid on
            handle_e(3) = plot(signal.t,e2,'b');
            handle_e(4) = plot(signal.t,e3,'r:');
            set(gca,'Xlim',[signal.t(1) signal.t(end)],'YScale','log','YLim',[0.01 10])
            set(gca,'TickLabelInterpreter','latex','fontsize',12);
            xlabel('Time [s]','interpreter','latex');
            ylabel('$$\Vert\mathbf{p} - \mathbf{\hat{p}}\Vert$$ [m]','interpreter','latex');
            p = get(gca,'position'); set(gca,'position', [p(1:3) 0.27])
            legend('RTI','EKF','RTSS','interpreter','latex','orientation','horizontal')
            
            PARAMS.VISUALIZE.handle_e = handle_e;
        end
        

        % plot CDF of RTI, EKF and RTSS errors
        if ~isempty(pos_idx)
            
            try
                handle_cdf = PARAMS.VISUALIZE.handle_cdf;
                [f_x,x] = ecdf(e1); set(handle_cdf(1),'Xdata',x,'Ydata',f_x);
                [f_x,x] = ecdf(e2); set(handle_cdf(2),'Xdata',x,'Ydata',f_x);
                [f_x,x] = ecdf(e3); set(handle_cdf(3),'Xdata',x,'Ydata',f_x);
            catch
                figure(3); clf; hold on; box on; grid on
                
                [f_x,x] = ecdf(e1); handle_cdf(1) = plot(x,f_x);
                [f_x,x] = ecdf(e2); handle_cdf(2) = plot(x,f_x);
                [f_x,x] = ecdf(e3); handle_cdf(3) = plot(x,f_x);

                set(gca,'Xlim',[.01 10],'XScale','log');
                set(gca,'TickLabelInterpreter','latex','fontsize',16);
                xlabel('$$\Vert\mathbf{p} - \mathbf{\hat{p}}\Vert$$ [m]','interpreter','latex');
                ylabel('F(x)','interpreter','latex');
                title('CDF','interpreter','latex')
                legend('RTI','EKF','RTSS','interpreter','latex')
                
                PARAMS.VISUALIZE.handle_cdf = handle_cdf;
            end

        end
        drawnow
    end
end