function ETASI_PlotResult(Result)
% Plot_ETASI_Result
% -------------------------------------------------------------------------
% Plot fitted results from Fit_ETASI_Hainzl.
%
% Input:
%   Result : output structure from Fit_ETASI_Hainzl
%
% Main plots:
%   1. True ETAS rate R0(t)
%   2. Apparent ETASI rate R(t)
%   3. Comparison between R0(t) and R(t)
%   4. Detection rate ratio R(t)/R0(t)
%   5. N0(t)=Tb*R0(t)
%   6. Magnitude-time catalog
%   7. Apparent ETASI magnitude density
%   8. Time-rescaling residual plot

if nargin < 1 || isempty(Result)
    error('Input Result structure is required.');
end

if ~isfield(Result, 'grid') || ~isfield(Result, 'event')
    error('Result must contain Result.grid and Result.event.');
end

grid = Result.grid;
event = Result.event;

% Use shifted time if available.
if isfield(Result, 'time0')
    eventTime = Result.time0;
    xLabelText = 'Time since first event';
else
    eventTime = Result.time;
    xLabelText = 'Time';
end

gridTime = grid.time;

% Magnitude for plotting
if isfield(Result, 'mag')
    magPlot = Result.mag;
elseif isfield(event, 'mag0')
    magPlot = event.mag0;
else
    magPlot = nan(size(eventTime));
end

% -------------------------------------------------------------------------
% Figure 1: true and apparent rates
% -------------------------------------------------------------------------
figure('Name', 'ETASI Rates', 'Color', 'w');
subplot(3,1,1);
plot(gridTime, grid.R0, 'r-', 'LineWidth', 1.8);
hold on;
% plot(eventTime, event.R0, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
ylabel('R_0(t)');
Fun_defaultAxes;
ax=gca; ax.YScale='log';

subplot(3,1,2);
plot(gridTime, grid.R, 'b-', 'LineWidth', 1.8);
hold on;
% plot(eventTime, event.R, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
ylabel('R(t)');
Fun_defaultAxes;
ax=gca; ax.YScale='log';

subplot(3,1,3);
plot(gridTime, grid.R0, 'r-', 'LineWidth', 1.8);
hold on;
plot(gridTime, grid.R, 'b-', 'LineWidth', 1.8);
xlabel(xLabelText);
ylabel('Rate');
Fun_defaultAxes;
ax=gca; ax.YScale='log';

% -------------------------------------------------------------------------
% Figure 2: detection ratio and N0
% -------------------------------------------------------------------------
figure('Name', 'ETASI Detection Effect', 'Color', 'w');
plot(gridTime, grid.detection_rate_ratio, 'k-', 'LineWidth', 2.0);
hold on;
% plot(eventTime, event.detection_rate_ratio, 'ro', 'MarkerSize', 3);
yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.2);
xlabel(xLabelText);
ylabel('R(t)/R_0(t)');
title('Detection rate ratio');
ylim([0, 1.05]);
Fun_defaultAxes;
% 
% subplot(2,1,2);
% plot(gridTime, grid.N0, 'm-', 'LineWidth', 1.5);
% hold on;
% plot(eventTime, event.N0, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
% 
% xlabel(xLabelText);
% ylabel('N_0(t) = T_b R_0(t)');
% title('Expected number of true events during blind time');
% legend('Grid N_0(t)', 'Event N_0(t_i)', 'Location', 'best');
% Fun_defaultAxes;

% -------------------------------------------------------------------------
% Figure 3: magnitude-time catalog with detection ratio
% -------------------------------------------------------------------------
figure('Name', 'Catalog and Detection Ratio', 'Color', 'w');
scatter(eventTime, magPlot, min(5.^(magPlot-min(magPlot))+12,200), event.detection_rate_ratio, 'filled');
xlabel(xLabelText);
ylabel('Magnitude');
cb = colorbar;
cb.Label.String = 'R(t)/R_0(t)';
colormap(flip(slanCM('viridis',30)));
Fun_defaultAxes;

% -------------------------------------------------------------------------
% Figure 5: time-rescaling residual check
% -------------------------------------------------------------------------
if isfield(Result, 'TransN') && isfield(Result, 'NumEvents')
    figure('Name', 'ETASI Time Rescaling', 'Color', 'w');
    tau = Result.TransN(:);
    numEvents = Result.NumEvents(:);
    plot(tau, numEvents, 'k', 'MarkerSize', 8, 'MarkerFaceColor', 'k','LineWidth',2);
    hold on;
    tauLine = linspace(min(tau), max(tau), 200)';
    plot(tauLine, tauLine, ':', 'LineWidth', 1.5,'Color',[.5 .5 .5]);
    xlabel('Transformed Time');
    ylabel('Cum.Events');
%     legend('Transformed events', 'Unit-rate Poisson expectation', ...
%         'Location', 'best');
    Fun_defaultAxes;
    axis('tight');
    set(gcf,'position',[300,50,700,700]);
end

end
