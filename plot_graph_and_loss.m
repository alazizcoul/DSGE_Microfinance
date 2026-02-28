%% ========================================================================
%  ANALYSE ET VISUALISATION GLOBALE - MODÈLE MICROFINANCE CÔTE D'IVOIRE
%  Auteur : AL AZIZ N'GOLO COULIBALY (Février 2026)
%% ========================================================================

clear ; close all; clc;

%% 1. Chargement des données,,
% -------------------------------------------------------------------------
if exist('results_comparison.mat', 'file')
    load('results_comparison.mat');
else
    error('Le fichier results_comparison.mat est introuvable. Lancez d''abord la simulation.');
end

%% 2. Paramétrage des graphiques
% -------------------------------------------------------------------------
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultLineLineWidth', 2);
T = T_periods;

%% 3. Variables à tracer
vars_as_points = {'zeta','vi','spr_b','RatioCap','ROA'};

vars_to_plot = {'zeta','vi','spr_b','B','K_m','RatioCap','ROA','output','consumption','investment'};
titles_vars_to_plot = {'a) Risque de crédit global','b) Norme de capitalisation','c) Marge d''intermédiation',...
                       'd) Prêts totaux','e) Fonds propres','f) Ratio de capitalisation','g) Rentabilité des actifs (ROA)',...
                       'h) Production ','i) Consommation','j) Investissement'};

%% 5. Tracer toutes les figures
%% 5. Tracer les figures principales

% --- Figure 1 : Choc de risque de crédit ---
plot_vars_4(vars_to_plot, titles_vars_to_plot, ...
            'NF_RC', 'ND_RC', ...
            'NF_RC_FP', 'ND_RC_FP', ...
            results, vars_as_points, ...
            'Figure I : Choc oositif de risque de crédit');

% --- Figure 2 : Choc de productivité ---
plot_vars(vars_to_plot, titles_vars_to_plot, ...
          'NF_Prod', 'ND_Prod', ...
          results, vars_as_points, ...
          'Figure II : Choc de productivité négative');

%% 6. Analyse numérique : fonction de perte sociale
fprintf('\n=============================================================\n');
fprintf('   ÉVALUATION DES PERFORMANCES (FONCTION DE PERTE SOCIALE)\n');
fprintf('=============================================================\n\n');

w_y  = 0; 
w_by = 1.0; 
w_vi = 1.0; 

scenarios_list = fieldnames(results);
fprintf('%-15s | %-12s | %-12s | %-10s\n', 'Scénario', 'Ecart-type(vi)', 'Ecart-type(B/Y)', 'Perte Tot.');
fprintf('-------------------------------------------------------------\n');

for s = 1:length(scenarios_list)
    scn = scenarios_list{s};
    if isfield(results.(scn), 'error'), continue; end

    irfs = results.(scn).irf_raw;
    L_y  = sum( (irfs.Y - results.(scn).ss_raw.Y).^2 ) / T;
    L_by  = sum( (irfs.B/irfs.Y - results.(scn).ss_raw.B/results.(scn).ss_raw.Y).^2 ) / T;
    L_vi = sum( (irfs.vi - results.(scn).ss_raw.vi).^2 ) / T;

    loss = (w_y * L_y) + (w_by * L_by) + (w_vi * L_vi);
    fprintf('%-15s | %-12.4f | %-12.4f | %-10.6f\n', scn, sqrt(L_vi), sqrt(L_by), loss);
end
fprintf('=============================================================\n');

%% 4. Fonction utilitaire pour tracer les variables
function plot_vars_4(vars, titles, scn1, scn2, scn3, scn4, results, vars_points, fig_title)

    figure('Name', fig_title, ...
           'Units','normalized', ...
           'Position',[0.1 0.1 0.9 0.8], ...
           'Color','w');

    n = length(vars);
    nrows = 2;
    ncols = ceil(n / nrows);

    T = length(results.(scn1).irf.(vars{1}));
    time = 1:T;

    annotation('textbox',[0 0.95 1 0.05], 'String', fig_title, ...
               'EdgeColor','none', ...
               'HorizontalAlignment','center', ...
               'FontWeight','bold', ...
               'FontSize',12);

    for i = 1:n
        v = vars{i};
        subplot(nrows, ncols, i); hold on; grid on;

        % Etat stationnaire
        if isfield(results.(scn1).ss, v)
            ss_val = results.(scn1).ss.(v);
            plot(time, ss_val * ones(1,T), 'k--', 'LineWidth',1);
        end

        % NF
        plot(time, results.(scn1).irf.(v), 'r-', 'LineWidth',1.0);
        % ND
        plot(time, results.(scn2).irf.(v), 'b-', 'LineWidth',1.0);
        % NF + FP
        plot(time, results.(scn3).irf.(v), 'r--', 'LineWidth',1.5);
        % ND + FP
        plot(time, results.(scn4).irf.(v), 'b--', 'LineWidth',1.5);

        title(titles{i}, 'FontWeight','bold');

        if ismember(v, vars_points)
            ylabel('Pct(%)');
        else
            ylabel('% Deviation ss');
        end

        if i == 1
            legend({'SS','NF','ND','NF+FP','ND+FP'}, 'Location','best');
        end
    end

    saveas(gcf, [regexprep(lower(fig_title),'[^a-z0-9]','_') '.png']);
end


function plot_vars(vars, titles, scn_NF, scn_ND, results, vars_points, fig_title)

    figure('Name', fig_title, ...
           'Units','normalized', ...
           'Position',[0.1 0.1 0.9 0.8], ...
           'Color','w');

    n = length(vars);
    nrows = 2;
    ncols = ceil(n / nrows);

    T = length(results.(scn_NF).irf.(vars{1}));
    time = 1:T;

    annotation('textbox',[0 0.95 1 0.05], 'String', fig_title, ...
               'EdgeColor','none', ...
               'HorizontalAlignment','center', ...
               'FontWeight','bold', ...
               'FontSize',12);

    for i = 1:n
        v = vars{i};
        subplot(nrows, ncols, i); hold on; grid on;

        % ===== ETAT STATIONNAIRE (NOIR) =====
        if isfield(results.(scn_NF).ss, v)
            ss_val = results.(scn_NF).ss.(v);
            plot(time, ss_val * ones(1,T), 'k--', 'LineWidth',1, ...
                 'DisplayName','Etat stationnaire');
        end

        % ===== NORME FIXE (ROUGE) =====
        if isfield(results.(scn_NF).irf, v)
            plot(time, results.(scn_NF).irf.(v), ...
                 'r*', 'LineWidth',0.5, ...
                 'DisplayName','Norme Fixe');
        end

        % ===== NORME DYNAMIQUE (BLEU) =====
        if isfield(results.(scn_ND).irf, v)
            plot(time, results.(scn_ND).irf.(v), ...
                 'b-', 'LineWidth',1, ...
                 'DisplayName','Norme Dynamique');
        end

        title(titles{i}, 'FontWeight','bold');

        if ismember(v, vars_points)
            ylabel('Pct(%)');
        else
            ylabel('% Deviation ss');
        end

        if i == 1
            legend('Location','best');
        end
    end

    % ===== Sauvegarde propre =====
    safe_title = lower(fig_title);
    safe_title = regexprep(safe_title, '[^a-z0-9]', '_');
    safe_title = regexprep(safe_title, '_+', '_');
    safe_title = regexprep(safe_title, '^_|_$', '');

    saveas(gcf, [safe_title '.png']);
end