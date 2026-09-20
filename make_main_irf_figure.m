clear
clc

%% Paths

folder = fileparts(mfilename('fullpath'));

results_folder = fullfile(folder,'results');
figures_folder = fullfile(folder,'figures');
%% Load results

% Baseline
S = load(fullfile(results_folder,'baseline.mat'),'results_repl');
baseline = S.results_repl;
clear S

% Pre-COVID
S = load(fullfile(results_folder,'ext_sample_1.precovid.mat'),'results_ext');
precovid = S.results_ext;
clear S

% Full sample with pandemic controls
S = load(fullfile(results_folder,...
    'ext_sample_3.pandemic_controls.mat'),'results_ext');
pandemic = S.results_ext;
clear S

%% Settings

hor = 24;
h = 0:hor;

row_names = { ...
    '1985Q1--2013Q2', ...
    '1985Q1--2019Q4', ...
    '1985Q1--2026Q2 + pandemic controls'};

column_names = { ...
    'GDP', ...
    'Year-ended inflation', ...
    'Cash rate', ...
    'Real TWI'};

units = {'%','ppt','ppt','%'};

%% Put models into cell arrays

models = {baseline, precovid, pandemic};

%% Common y-axis ranges by variable
% Important: same scale down each column makes comparison meaningful.

ylims = NaN(4,2);

for v = 1:4

    allvals = [];

    for r = 1:3

        M = models{r};

        switch v
            case 1
                point = M.irf_gdp(:);
                ci = M.gdp_CI;

            case 2
                point = M.irf_inf(:);
                ci = M.inf_CI;

            case 3
                point = M.irf_rate(:);
                ci = M.rate_CI;

            case 4
                point = M.irf_rtwi(:);
                ci = M.rtwi_CI;
        end

        % Point estimate and outer 95% bands
        allvals = [allvals; point; ci(:,1); ci(:,4)];

    end

    ymin = min(allvals);
    ymax = max(allvals);

    pad = 0.08*(ymax-ymin);

    if pad == 0
        pad = 0.1;
    end

    ylims(v,:) = [ymin-pad ymax+pad];

end

%% Plot

fig = figure( ...
    'Units','centimeters', ...
    'Position',[2 2 19 18], ...
    'Color','w');

t = tiledlayout(3,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for r = 1:3

    M = models{r};

    point_irfs = { ...
        M.irf_gdp, ...
        M.irf_inf, ...
        M.irf_rate, ...
        M.irf_rtwi};

    cis = { ...
        M.gdp_CI, ...
        M.inf_CI, ...
        M.rate_CI, ...
        M.rtwi_CI};

    for v = 1:4

        ax = nexttile;

        point = point_irfs{v};
        ci = cis{v};

        hold(ax,'on')

        % 95% confidence intervals
        plot(ax,h,ci(:,1),'--r','LineWidth',0.8)
        plot(ax,h,ci(:,4),'--r','LineWidth',0.8)

        % 85% confidence intervals
        plot(ax,h,ci(:,2),'--b','LineWidth',0.8)
        plot(ax,h,ci(:,3),'--b','LineWidth',0.8)

        % Point estimate
        plot(ax,h,point,'-k','LineWidth',1.5)

        % Zero line
        yline(ax,0,'-k','LineWidth',0.5)

        xlim(ax,[0 hor])
        ylim(ax,ylims(v,:))

        box(ax,'off')
        grid(ax,'off')

        ax.FontSize = 8;
        ax.TickDir = 'out';

        % Variable headings only along first row
        if r == 1
            title(ax,column_names{v}, ...
                'FontSize',9, ...
                'FontWeight','normal')
        end

        % X labels only on bottom row
        if r == 3
            xlabel(ax,'Quarters','FontSize',8)
        else
            ax.XTickLabel = [];
        end

        % Unit labels
        ylabel(ax,units{v},'FontSize',8)

        % Row label on first panel of each row
        if v == 1
            text(ax,-0.29,0.5,row_names{r}, ...
                'Units','normalized', ...
                'Rotation',90, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',8);
        end

    end
end
%% Overall title

title(t, ...
    'Impulse Responses to a 10 percent exogenous real exchange rate appreciation', ...
    'FontSize',10, ...
    'FontWeight','bold')

%% Export

exportgraphics(fig, ...
    fullfile(figures_folder,'main_irfs_3x4.pdf'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(figures_folder,'main_irfs_3x4.png'), ...
    'Resolution',300);