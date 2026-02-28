%% ========================================================================
%  COMPARAISON DES NORMES DE CAPITALISATION : FIXE vs DYNAMIQUE
%  Auteur : AL AZIZ N'GOLO COULIBALY (Adaptation Microfinance)
%  Mise à jour : Correction EXP() et Choc Permanent e_vi
%% ========================================================================

clear ; close all; clc;

%% -------------------------------------------------------------------------
%  1. CONFIGURATION DES SCÉNARIOS
%% -------------------------------------------------------------------------
scenarios = {
    'NF_RC',    'static',  {'e_zeta_e', 'e_zeta_i', 'e_vi'},        'Norme Fixe (Choc de Norme) + Choc Risque Crédit';
    'ND_RC',    'dynamic', {'e_zeta_e', 'e_zeta_i'},                'Norme Dynamique + Choc Risque Crédit';
    'NF_Prod',  'static',  {'e_A_e', 'e_vi'},                       'Norme Fixe (Choc de Norme) + Choc Productivité';
    'ND_Prod',  'dynamic', {'e_A_e'},                               'Norme Dynamique + Choc Productivité';
    'NF_RC_FP', 'static',  {'e_zeta_e', 'e_zeta_i', 'e_vi', 'e_eps_K_m'}, 'Norme Fixe (Choc de Norme) + Choc RC + Injection FP';
    'ND_RC_FP', 'dynamic', {'e_zeta_e', 'e_zeta_i', 'e_eps_K_m'},   'Norme Dynamique + Choc RC + Injection FP';
   % 'NF_RC_FP2', 'static',  {'e_zeta_e', 'e_zeta_i', 'e_vi', 'e_eps_K_m'}, 'Norme Fixe (Choc de Norme) + Choc RC + Injection FP';
   %  'ND_RC_FP2', 'dynamic', {'e_zeta_e', 'e_zeta_i'},   'Norme Dynamique + Choc RC Vs RC + Injection FP';
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
};


%% -------------------------------------------------------------------------
%  2. MAGNITUDE DES CHOCS
%% -------------------------------------------------------------------------
shock_magnitude = struct();
shock_magnitude.e_vi       =  1.0;  % Choc permanent (ex: +0.25 point)
shock_magnitude.e_zeta_e   =  1.0;   
shock_magnitude.e_zeta_i   =  1.0;   
shock_magnitude.e_A_e      = -1.0;   
shock_magnitude.e_eps_K_m  = -1.0;   

T_periods = 40;

%% -------------------------------------------------------------------------
%  4 & 5. VARIABLES ET INITIALISATION
%% -------------------------------------------------------------------------
%irf_vars = {'Y', 'C', 'I', 'B', 'K_m', 'vi', 'r_be', 'r_bh', 'zeta_e', 'zeta_i', 'pie'};

irf_vars = {'Y', 'C', 'I', 'B', 'K_m', 'vi', 'r_be', 'r_bh', 'zeta_e', 'zeta_i', 'zeta', 'pie','RatioCap', 'spr_b' ...
            'output', 'consumption', 'investment', 'loansH', 'loansF', 'deposits',...
            'IMFcapital', 'interestPol', 'interestH', 'interestF', 'interestDep', 'inflation', 'ROA'};
% Variables à traiter linéairement (déjà en niveau ou taux dans le .mod)
vars_no_exp = {'RatioCap','ROA', 'interestPol', 'interestH', 'interestF', 'interestDep', ...
               'inflation', 'output', 'consumption', 'investment', 'loansH', ...
               'loansF', 'deposits', 'IMFcapital', 'pie'};

% Variables à afficher en points de pourcentage (multipliées par 100)
vars_as_points = {'vi', 'r_be', 'r_bh','spr_b',  'zeta','pie', 'ROA', 'interestPol', ...
                  'interestH', 'interestF', 'interestDep', 'inflation','RatioCap'};

% Variables à afficher en pourcentage de déviation par rapport à l'état stationnaire)
vars_as_percents = {'Y', 'C', 'I', 'B', 'K_m','output', 'consumption', 'investment', ...
                    'loansH', 'loansF', 'deposits', 'IMFcapital', };
              
              
results = struct();

%% -------------------------------------------------------------------------
%  6. BOUCLE PRINCIPALE
%% -------------------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('  SIMULATION DES SCÉNARIOS (MODE MICROFINANCE)\n');
fprintf('==================================================\n\n');

for s = 1:size(scenarios, 1)
    scenario_name = scenarios{s, 1};
    norm_type     = scenarios{s, 2};
    shocks_list   = scenarios{s, 3};
    description   = scenarios{s, 4};

    fprintf('Scénario %d/%d : %s\n', s, size(scenarios, 1), description);

    try
        %% instruction : Ajouter un code qui vérifie si M_ , oo_ et options_ et les supprimer avant le eval
        eval(sprintf('dynare model_microf -Dnorm="%s" noclearall', norm_type));
        M_saved = M_; oo_saved = oo_; options_saved = options_;
        results.(scenario_name).norm = norm_type;
        ss_vec = oo_saved.dr.ys;

        %% -- 6c. Construction de la matrice de chocs --
        ex_ = zeros(T_periods, M_saved.exo_nbr);
        for shock_index = 1:length(shocks_list)
            shock_name = shocks_list{shock_index};
            exo_index = find(strcmp(cellstr(M_saved.exo_names), shock_name));
            if isempty(exo_index), continue; end

            shock_std = sqrt(M_saved.Sigma_e(exo_index, exo_index));
            val = shock_magnitude.(shock_name) * shock_std;

            if strcmp(shock_name, 'e_vi')
                ex_(:, exo_index) = val; % Permanent (norme statique à 15% alors que le niveau stat est de 12.5% )
            elseif strcmp(shock_name, 'e_eps_K_m')
                ex_(2, exo_index) = val; % Choc de fnds propre à la deuxième période (réaction au choc de risque)               
            else
                ex_(1, exo_index) = val; % les autres chocs y compris e_A_e et e_zeta(e/i) sont à la période 1
            end
        end

        %% -- 6d. Simulation --
        y_ = simult_(M_saved, options_saved, ss_vec, oo_saved.dr, ex_, 1);

        %% -- 6e. Sauvegarde simplifiée des IRF --
        for v = 1:length(irf_vars)

            vname = irf_vars{v};
            idx = find(strcmp(cellstr(M_saved.endo_names), vname));
            if isempty(idx), continue; end

            ss_raw   = ss_vec(idx);
            path_raw = y_(idx, 2:end);

            % =====================================================
            % 1. Déterminer si exp() est nécessaire
            % =====================================================
            if ismember(vname, vars_no_exp)
                ss_val   = ss_raw;
                path_val = path_raw;
            else
                ss_val   = exp(ss_raw);
                path_val = exp(path_raw);
            end

            % =====================================================
            % 2. Calcul de la déviation
            % =====================================================
            deviation = path_val - ss_val;

            % =====================================================
            % 3. Règle d’affichage
            % =====================================================
            % =====================================================
            % 3. Règle d'affichage
            % =====================================================
            if ismember(vname, vars_as_points)

                results.(scenario_name).irf.(vname)       = path_val * 100;
                results.(scenario_name).ss.(vname)        = ss_val * 100;
                % VERSION NIVEAU BRUT
                results.(scenario_name).irf_raw.(vname)   = path_val;
                results.(scenario_name).ss_raw.(vname)    = ss_val;

            else

                results.(scenario_name).irf.(vname)       = (deviation ./ ss_val) * 100;
                results.(scenario_name).ss.(vname)        = 0;
                % VERSION NIVEAU BRUT
                results.(scenario_name).irf_raw.(vname)   = path_val;
                results.(scenario_name).ss_raw.(vname)    = ss_val;

            end
            
        end    

        fprintf('  [OK] Simulation réussie\n\n');
    catch ME
        fprintf('  [ERREUR] %s\n\n', ME.message);
    end
end 

save('results_comparison.mat', 'results', 'scenarios', 'shock_magnitude', 'T_periods');
fprintf('Calculs terminés. Les déviations utilisent désormais exp() pour les variables concernées.\n');