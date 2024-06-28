%% Analysis code for: The status of vernier acuity following late sight 
%  onset (Vogelsang et al., under review) 

%% PART 1: SETUP
%% Folders and paths
clear all; close all; 
results_dir = 'figures';
mkdir(results_dir); 

%% Extract data
load("datafiles.mat");
n_prakash_timepoints = size(prakash_vernier_all, 2); % Longitudinal timepoints
n_prakash = size(prakash_vernier_all, 1); % Nr. of Prakash patients

% Prakash data: mean + standard error
prakash_resacuity_mean = mean(prakash_resacuity_all, 'omitnan');
prakash_vernier_mean = mean(prakash_vernier_all, 'omitnan');
prakash_resacuity_se = std(prakash_resacuity_all, 'omitnan')./sqrt(sum(~isnan(prakash_resacuity_all)));
prakash_vernier_se = std(prakash_vernier_all, 'omitnan')./sqrt(sum(~isnan(prakash_vernier_all)));

% Control data: mean + standard error
control_200_vernier_mean = mean(control_200_vernier, 'omitnan');
control_200_vernier_se = std(control_200_vernier, 'omitnan')./sqrt(sum(~isnan(control_200_vernier)));
control_500_vernier_mean = mean(control_500_vernier, 'omitnan');
control_500_vernier_se = std(control_500_vernier, 'omitnan')./sqrt(sum(~isnan(control_500_vernier)));

%% PART 2: MAKE FIGURES
% Constants for figures
figFontSize = 12;
plot_colors_vec = [70, 6, 23; 
                  157, 2, 8;
                  210, 10, 2;
                  226, 60, 4;
                  242, 130, 6;
                  255 186 8;
                  70, 150, 230]./255;
plot_gray_vec = [0.4, 0.4, 0.4];

%% Figure 1
fig = figure('Position', [20, 70, 800, 600]); 

% Panel C (Resolution Acuity) & Panel D (Vernier Acuity)
acuity_mean = {prakash_resacuity_mean, prakash_vernier_mean};
acuity_se = {prakash_resacuity_se, prakash_vernier_se};
acuity_names = {'Resolution', 'Vernier'};
acuity_units = {'(c/deg)', '(1/min)'};
for i = 1:2
    subplot(2, 2, i);
    b = bar(1:n_prakash_timepoints, acuity_mean{i}, 'Linewidth', 1.5);
    b.FaceColor = 'flat';
    b.CData = plot_colors_vec(1:n_prakash_timepoints, :);
    hold on; 
    errorbar(1:n_prakash_timepoints, acuity_mean{i}, acuity_se{i}, ...
            'x', 'Color', plot_gray_vec, 'Linewidth', 1.5);
    title(strcat(acuity_names{i}, ' acuity means'));
    ylabel(strcat(acuity_names{i}, " acuity ", acuity_units{i}));
    xlabel('Timepoint')
    ax = gca; 
    ax.FontSize = figFontSize;
end

% Panel E: Scatter plot
subplot(2, 2, 3); 
scatter(prakash_resacuity_mean, prakash_vernier_mean, 80, ...
    plot_colors_vec(1:n_prakash_timepoints, :), 'filled', 'Linewidth', 2);
xlabel('Resolution acuity (c/deg)');
ylabel('Vernier acuity (1/min)'); 
title(strcat("Longitudinal means"))
ax = gca; 
ax.FontSize = figFontSize; 
ax.XLim(1) = 0; 
ax.YLim(1) = 0;
grid on;
hold on;
plot([0, 15], [0, 0.5], '--', 'Color', plot_gray_vec, 'LineWidth', 2);
xlim([0 15])

% Panel F: Individual data at timepoint 6
subplot(2, 2, 4);
t = 6;
scatter(prakash_resacuity_all(:, t), prakash_vernier_all(:, t), 40, ...
    plot_colors_vec(t, :), 'Linewidth', 2);
title("Indiv. data at final timepoint");
xlabel('Resolution acuity (c/deg)');
ylabel('Vernier acuity (1/min)');
ax = gca;
ax.FontSize = figFontSize;
ax.XLim = [0, ceil(max(prakash_resacuity_all(:)))]; 
ax.YLim = [0, ceil(max(prakash_vernier_all(:)))];
grid on;
hold on;
plot([0, 30], [0, 1], '--', 'Color', plot_gray_vec, 'LineWidth', 2);
xlim([0 30]);

saveas(fig, fullfile(results_dir, 'figure1.png'));

%% Supplementary Figure 1
fig = figure('Position', [20, 70, 600, 700]);
colorResacuity = [0.35 0.35 0.35]; 
colorVernier = [63, 152, 43]./255;

% For each participant:
for i = 1:n_prakash
    subplot(5, 2, i);
    
    % First y-axis (left side)
    yyaxis left;
    plot(1:n_prakash_timepoints, prakash_resacuity_all(i, :), 'o-', 'Color', ...
        colorResacuity, 'LineWidth', 2, 'MarkerFaceColor', colorResacuity);
    ylabel('Resolution');
    ylim([0 5]);
    yticks([0 5]);
    set(gca, 'ycolor', colorResacuity); 
    
    % Second y-axis (right side)
    yyaxis right;
    plot(1:n_prakash_timepoints, prakash_vernier_all(i, :), 'd-', 'Color', ...
        colorVernier, 'LineWidth', 2, 'MarkerFaceColor', colorVernier);
    ylabel('Vernier');
    ylim([0 1]);
    set(gca, 'ycolor', colorVernier); 
    yticks([0 1]);

    % Handling NaNs: plotting dashed lines between points adjacent to NaNs
    hold on;
    for k = 1:n_prakash_timepoints-1
        if isnan(prakash_resacuity_all(i, k)) && ~isnan(prakash_resacuity_all(i, k+1))
            start_idx = k;
            while isnan(prakash_resacuity_all(i, start_idx)) && start_idx > 1
                start_idx = start_idx - 1;
            end
            if start_idx < k 
                yyaxis left;
                plot([start_idx k+1], prakash_resacuity_all(i, [start_idx k+1]), ...
                    '--', 'Color', colorResacuity, 'LineWidth', 2);
            end
        end
        if isnan(prakash_vernier_all(i, k)) && ~isnan(prakash_vernier_all(i, k+1))
            start_idx = k;
            while isnan(prakash_vernier_all(i, start_idx)) && start_idx > 1
                start_idx = start_idx - 1;
            end
            if start_idx < k 
                yyaxis right;
                plot([start_idx k+1], prakash_vernier_all(i, [start_idx k+1]), ...
                    '--', 'Color', colorVernier, 'LineWidth', 2);
            end
        end
    end
    xlim([1 n_prakash_timepoints]);
    xticks(1:n_prakash_timepoints);
    if i == 9 || i == 10
        xlabel('Timepoints');
    else
        xlabel('');
    end
title(strcat("Participant ",num2str(i)))
ax = gca; ax.FontSize = 11;
end

saveas(fig, fullfile(results_dir, 'supplementary_fig1.png'));

%% Supplementary Figure 2
fig = figure('Position', [20, 70, 1250, 700]); 
for t = 1:n_prakash_timepoints
    corrval = corrcoef(prakash_vernier_all(:, t), ...
        prakash_resacuity_all(:, t), 'rows', 'complete');
    subplot(2, 3, t);
    scatter(prakash_resacuity_all(:, t), prakash_vernier_all(:, t), ...
        40, plot_colors_vec(t, :), 'Linewidth', 2);
    hold on;
    title(strcat("Timepoint ", num2str(t), " (r = ", num2str(corrval(2)),")"))
    xlabel('Resolution acuity (c/deg)');
    ylabel('Vernier acuity (1/min)');
    ax = gca;
    ax.FontSize = figFontSize+1;
    ax.XLim = [0, ceil(max(max(prakash_resacuity_all)))]; 
    ax.YLim = [0, ceil(max(max(prakash_vernier_all)))];
    grid on;
    plot([ax.XLim(1), ax.XLim(2)], [ax.XLim(1) ./ 30, ax.XLim(2) ./ 30], ...
        '--', 'Color', plot_gray_vec, 'LineWidth', 2);
end
saveas(fig, fullfile(results_dir, 'supplementary_fig2.png'));

%% Figure 2
fig = figure('Position', [20, 70, 900, 600]); 

% Panel A: Plot 20/500 comparison
subplot(2, 3, 1);
t = 3;
b = bar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1,:) = plot_colors_vec(t, :); 

errorbar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], ...
    [prakash_vernier_se(t), control_500_vernier_se], 'x', 'Color', plot_gray_vec, 'Linewidth', 1.5)
title('1.2 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 1];
ax.XTickLabel={'Prakash','Control'};

% Panel C: Plot 20/200 comparison
subplot(2, 3, 4);
t = 6;
b = bar(1:2, [prakash_vernier_mean(t), control_200_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1, :) = plot_colors_vec(t, :);
b.CData(2, :) = plot_colors_vec(t + 1, :); 

errorbar(1:2, [prakash_vernier_mean(t), control_200_vernier_mean], ...
    [prakash_vernier_se(t), control_200_vernier_se], 'x', 'Color', plot_gray_vec, 'Linewidth', 1.5)
title('3 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 1];
ax.XTickLabel = {'Prakash','Control'};

% Panel B: Individual data (20/500 controls + Prakash at timepoint 3)
subplot(2, 3, 2:3);
t = 3;
non_nan_controldata = control_500_vernier(~isnan(control_500_vernier));
control_length = length(non_nan_controldata);
prakash_non_nan_idx = ~isnan(prakash_vernier_all(:, t));
prakash_length = sum(prakash_non_nan_idx);
b = bar([prakash_vernier_all(prakash_non_nan_idx, t); 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(plot_colors_vec(t, :), prakash_length, 1); 
ax = gca; 
ax.FontSize = figFontSize; 
ax.YLim = [0, 1];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

% Panel B: Individual data (20/200 controls + Prakash at timepoint 6)
subplot(2, 3, 5:6);
t = 6;
non_nan_controldata = control_200_vernier(~isnan(control_200_vernier));
control_length = length(non_nan_controldata);
prakash_non_nan_idx = ~isnan(prakash_vernier_all(:, t));
prakash_length = sum(prakash_non_nan_idx);
b = bar([prakash_vernier_all(prakash_non_nan_idx, t); 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(plot_colors_vec(t, :), prakash_length, 1); 
b.CData(prakash_length + 2:end, :) = repmat(plot_colors_vec(t + 1, :), control_length, 1);
ax = gca; 
ax.FontSize = figFontSize;
ax.YLim = [0, 1];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

saveas(fig, fullfile(results_dir, 'figure2.png'));

%% Supplementary Figure 3
fig = figure('Position', [20, 70, 900, 300]); 

% Panel A: Plot 20/500 comparison (but with Prakash timepoint 2, not 3)
subplot(1, 3, 1);
t = 2;
b = bar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], 'Linewidth', 1.5);
hold on; 
b.FaceColor = 'flat';
b.CData(1,:) = plot_colors_vec(t, :); 

errorbar(1:2, [prakash_vernier_mean(t), control_500_vernier_mean], ...
    [prakash_vernier_se(t), control_500_vernier_se], 'x', 'Color', plot_gray_vec, 'Linewidth', 1.5)
title('1.2 c/deg comparison');
ylabel('Vernier acuity (1/min)'); 
ax = gca;
ax.FontSize = figFontSize;
ax.YLim = [0, 1];
ax.XTickLabel={'Prakash','Control'};

% Panel B: Individual data (20/500 controls + Prakash at timepoint 2)
subplot(1, 3, 2:3);
t = 2;
non_nan_controldata = control_500_vernier(~isnan(control_500_vernier));
control_length = length(non_nan_controldata);
prakash_non_nan_idx = ~isnan(prakash_vernier_all(:, t));
prakash_length = sum(prakash_non_nan_idx);
b = bar([prakash_vernier_all(prakash_non_nan_idx, t); 0; non_nan_controldata']);
b.FaceColor = 'flat'; 
b.CData(1:prakash_length, :) = repmat(plot_colors_vec(t, :), prakash_length, 1); 
ax = gca; 
ax.FontSize = figFontSize; 
ax.YLim = [0, 1];
ylabel('Vernier acuity (1/min)'); 
title('Individual participants');

saveas(fig, fullfile(results_dir, 'supplementary_fig3.png'));

%% Figures 3A&B
% Coordinates for axes
x_axis = [8 10 12 14 16 18 20 22];
y_axis_left = [1.25 2.5 5 10 20 40];
y_axis_right = [12 6 3 1.5 0.75];

% Figure 3A: Extracted data from Shimojo & Held:
left_x = [9.4555   11.4853   13.5397   15.5333   17.6876   19.6402];
left_y = [ 46.1868   28.4896   16.9638   10.3286    5.3799    6.1500];
right_x = [9.5006   11.4931   13.5213   15.4595   17.6036   19.6269];
right_y = [2.0501    1.9825    2.5257    2.3707    2.7363    3.2117];

% Figure 3A: Left y-axis
fig = figure('Position',[100 100 450 400]);
scatter(log(left_x), log(left_y), 120, 'd', 'MarkerEdgeColor', ...
    [0, 0, 0]./255, 'MarkerFaceColor', [63, 152, 43]./255);
hold on;
set(gca, 'YDir', 'reverse'); % Invert axis
scatter(log(right_x), log(30./right_y), 120, 'MarkerEdgeColor', ...
    [0, 0, 0]./255, 'MarkerFaceColor', [80, 53, 162]./255);
x = 8:2:22;
xticks(log(x));
xticklabels(arrayfun(@num2str, x, 'UniformOutput', false));
xlim([log(8), log(22)]);
y = y_axis_left;
yticks(log(y));
yticklabels(arrayfun(@num2str, y, 'UniformOutput', false));
ylim([log(1.25), log(70)])

% Figure 3A: right y-axis
yyaxis right
set(gca, 'YDir', 'reverse'); % invert
ax = gca;
ax.YColor = 'k';
y = y_axis_left;
yticks(log(y));
yticklabels(arrayfun(@num2str, 30./y, 'UniformOutput', false));
ylim([log(1.25), log(70)])
legend({'Vernier acuity','Resolution acuity'}, 'FontSize',12,'Location','southeast');

% Figure 3B: extract data from Prakash data
x = [1 4 53 113 330 1530];
left_y = 1./prakash_vernier_mean;
right_y = prakash_resacuity_mean;

% Figure 3B: left y-axis
yyaxis left
fig2 = figure('Position',[100 100 450 400]);
scatter(log(x), log(left_y), 120, 'd', 'MarkerEdgeColor', ...
    [0, 0, 0]./255, 'MarkerFaceColor', [63, 152, 43]./255);
set(gca, 'YDir', 'reverse'); % Invert y-axis
hold on;
scatter(log(x), log(30./right_y), 120, 'MarkerEdgeColor', ...
    [0, 0, 0]./255, 'MarkerFaceColor', [80, 53, 162]./255);
xticks(log(x))
xticklabels(arrayfun(@num2str, x, 'UniformOutput', false));
y = y_axis_left;
yticks(log(y));
yticklabels(arrayfun(@num2str, y, 'UniformOutput', false));
ylim([log(1.25), log(70)])

% Figure 3B: Create secondary y-axis
yyaxis right
set(gca, 'YDir', 'reverse'); % Invert y-axis
ax = gca;
ax.YColor = 'k';
y = y_axis_left;
yticks(log(y));
yticklabels(arrayfun(@num2str, 30./y, 'UniformOutput', false));
ylim([log(1.25), log(70)])
legend({'Vernier acuity','Resolution acuity'}, 'FontSize',12,'Location','southeast');
yyaxis left

% Save both figures
saveas(fig, fullfile(results_dir, 'figure3a.png'));
saveas(fig2, fullfile(results_dir, 'figure3b.png'));

%% PART 3: STATISTICS
% Stats for group comparison (Figure 2): 
disp('Stats for group comparisons (Figure 2):');
[a, b, c, d] = ttest2(prakash_vernier_all(:, 2), control_500_vernier',"Tail","left"); 
disp(strcat('Prakash timepoint 2 vs. Control 20/500: t(', ...
    num2str(d.df),') = ', num2str(d.tstat), ', p = ', num2str(b)));
[a, b, c, d] = ttest2(prakash_vernier_all(:, 3), control_500_vernier',"Tail","left"); 
disp(strcat('Prakash timepoint 3 vs. Control 20/500: t(', ...
    num2str(d.df),') = ', num2str(d.tstat), ', p = ', num2str(b)));
[a, b, c, d] = ttest2(prakash_vernier_all(:, 6), control_200_vernier',"Tail","left"); 
disp(strcat('Prakash timepoint 6 vs. Control 20/200: t(', ...
    num2str(d.df),') = ', num2str(d.tstat), ', p = ', num2str(b)));
disp(' ');

% Stats for the line comparisons (Figure 1E): 
disp('Stats for line comparisons (Figure 1E):');
for timept = 1:n_prakash_timepoints
    [a,b,c,d] = ttest(prakash_vernier_all(:,timept), prakash_resacuity_mean(timept)./30, "Tail", "right");
    disp(strcat('Vernier acuity > equivalence line, for timepoint ', ...
        num2str(timept), ': t(', num2str(d.df),') = ', ...
        num2str(d.tstat), ', p = ', num2str(b), ', p(corr) = ', num2str(b*6) ));
end
disp(' ');

% Stats for correlation (Figure 1E): 
disp('Stats for correlation (Figure 1E):');
res_ver_correlation = corrcoef(prakash_resacuity_mean, prakash_vernier_mean);
disp(strcat('Correlation between the two acuity means over time: ', ...
    num2str(res_ver_correlation(2))));