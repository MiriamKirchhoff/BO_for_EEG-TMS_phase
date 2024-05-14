%% Simulation main script
% 1D search
% Code by Miriam Kirchhoff, main script and functions adapted from
% (1) Aino Nieminen & Jukka Sarvas, Aalto University
% (2) Martin Krasser (https://github.com/krasserm)
% v08072022
%
% version   1.0, 19.04.2024
% author    Miriam Kirchhoff
% project   C2B

clear
close all

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings = BO_settings;
settings.type_acquisition = {'probEI', 'EI', 'random'};
tic

%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load or generate data

% if path was specified, load data
if settings.path_load
    disp('Loading data')
    load(settings.path_load)
    disp('Loading data completed')
else
    % if path was not specified, generate random data
    disp('Generate example data')
    [all_results, all_trials] = generate_example_data(settings);
    disp('Generating data completed')
end

current_datetime = strrep(datestr(now),':','_');


%% Generate computation grid

% number of potential stimulation targets
init.T_N = 200;

% min and max of stimulation targets
init.T_min = -pi;
init.T_max = pi;

% number of targets in the grid search
grid_number = 8;
T = linspace(init.T_min, init.T_max, init.T_N);
n = numel(T);
iter_max = settings.n_iterations + 1;


%% Initialize result arrays

% store dimensions of results arrays
dim = ["reps", "iter", "acq", "mod", "sub"];
simulation_results_opt_estimate = NaN(settings.n_repetitions, settings.n_iterations+ 1, length(settings.type_acquisition), length(settings.type_model), length(settings.subjects));
simulation_results_error        = simulation_results_opt_estimate;


%% Initialize figure

if settings.plotting
    online_plot = figure('Name','Measured Data', 'NumberTitle','off');
end % if settings.plotting

%% Calculate ground truth for simulation

% iterate over subject datasets
for subject_idx = 1:length(settings.subjects)
    current_subject = settings.subjects(subject_idx);

    % generate ground truth for computation grid using sliding window
    % method

    % settings
    subject_data.phases = all_trials.phase(...
        all_trials.participant == current_subject & ~all_trials.outlier_all);
    subject_data.transformed_mep_amplitudes = all_trials.MEP_log(...
        all_trials.participant == current_subject & ~all_trials.outlier_all)';
    simulation_parameters.percentage_data = 0.1;          % Percentage of complete window to be considered
    simulation_parameters.steps =   n;
    simulation_parameters.wintype = 'winsize';
    simulation_parameters.winsize = (init.T_max - init.T_min)*simulation_parameters.percentage_data;  % window size of simulation
    simulation_parameters.limits =  [init.T_min, init.T_max];

    % calculation
    ground_truth = simulation_ground_truth_slidewindow(simulation_parameters, subject_data, 0, false, false, false);
    clear simulation_parameters subject_data


    %% Start loops

    % iterate over selected models
    for model_idx = 1:length(settings.type_model)
        current_model = settings.type_model{model_idx};

        % iterate over selected acquisition functions
        for acquisition_idx = 1:length(settings.type_acquisition)
            current_acquisition = settings.type_acquisition{acquisition_idx};

            % this for can be replaced by parallel processing using parfor
            % if plotting == 0
            parfor repetition_current = 1:settings.n_repetitions

                % output message: current repetition
                fprintf(['start repetition %g/%g, acquisition %g/%g, ' ...
                    'model %g/%g, subject %g/%g \n'], ...
                    repetition_current, settings.n_repetitions, ...
                    acquisition_idx, length(settings.type_acquisition),...
                    model_idx, length(settings.type_model), ...
                    subject_idx, length(settings.subjects))


                %% DATA GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %% Initialize

                t = [];
                X_ind = [];
                acq = [];

                % Get initial measurement locations
                temp = floor(length(T)/settings.n_init_samples);

                % select random point in first 1/n area
                X_new = randi(temp);

                % add equally spaced samples
                for i = 1:settings.n_init_samples-1
                    X_new = [X_new X_new(1) + i*temp];
                end % for i = 1:settings.n_init_samples-1
                % clear i temp

                %% SIMULATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for iteration = 1:iter_max

                    %% Get simulated target values t

                    t_new = [];
                    % sample all required new datapoints
                    for i = 1:length(X_new)
                        t_new = [t_new normrnd(ground_truth.mean(X_new(i)), ...
                            ground_truth.sd(X_new(i)))];
                    end % for i = 1:length(X_new)

                    % Append to data arrays
                    t = [t t_new];
                    X_ind = [X_ind X_new];
                    X = T(X_ind);
                    % clear t_new X_new


                    %% FIT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % select which model to calculate
                    switch current_model
                        case 'GP'
                            if iteration == 1
                                GP_fit = struct('postmean', [], 'postvar', [], ...
                                    'postmeanX', [], 'kernel', settings.type_kernel);
                            end
                            [y, y_var, GP_fit] = BO_fit_GP(X_ind, t, T, GP_fit);

                        case 'bayesian regression'
                            % fit a bayesian linear regression
                            [m_N, S_N, beta] = blr_regression(X', t');
                            phi_test = [ones(size(T')), sin(T'), cos(T')];
                            [y, y_var] = blr_posterior_predictive(phi_test, m_N, S_N, beta);
                    end % switch current_model


                    %% ACQUISITION FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Select which acquisition function to fit
                    switch current_acquisition
                        case 'KG'
                            % Choose the next stimulation parameter with the knowledge gradient.
                            % Works well in 1D, but in higher dimensional cases,
                            % computationally cheaper acquisition function may be needed.
                            switch current_model
                                case 'GP'
                                    [X_new,acq] = tms_nextpoint_kg(X,t',GP_fit.mu0,...
                                        GP_fit.a0,GP_fit.a1,GP_fit.lam2,T,GP_fit.muz3,...
                                        ones(n,1)*GP_fit.lam2(1),GP_fit.kernel);

                                otherwise
                                    [X_new,acq] = blr_nextpoint_kg(X', t', beta, m_N, S_N, T');

                                    % if acquisition function is noncontinuous
                                    if sum(isnan(acq))
                                        X_new = randi(length(T));
                                        acq = NaN(size(T));
                                        disp('selecting random simulation target due to nonconvergence')
                                    end
                            end

                        case 'probKG'
                            % Choose the next stimulation parameter with the knowledge gradient.
                            % Works well in 1D, but in higher dimensional cases,
                            % computationally cheaper acquisition function may be needed.
                            switch current_model
                                case 'GP'
                                    [X_new,acq] = tms_nextpoint_kg(X,t',GP_fit.mu0,...
                                        GP_fit.a0,GP_fit.a1,GP_fit.lam2,T,GP_fit.muz3,...
                                        ones(n,1)*GP_fit.lam2(1),GP_fit.kernel);

                                otherwise
                                    [X_new,acq] = blr_nextpoint_kg(X', t', beta, m_N, S_N, T');

                                    % if acquisition function is noncontinuous
                                    if sum(isnan(acq))
                                        X_new = randi(length(T));
                                        acq = NaN(size(T));
                                        disp('selecting random simulation target due to nonconvergence')
                                    end

                                    acq = acq';
                            end

                            % if acquisition function is noncontinuous
                            if sum(isnan(acq))
                                X_new = randi(length(T));
                                acq = NaN(size(T));
                                disp('selecting random simulation target due to nonconvergence')
                            else

                                % select the next target randomly from the acq,
                                % interpreting the acq as the probability of
                                % the target being selected

                                Y_new = randsrc(1,1,[T ; (acq/sum(acq))']);
                                X_new = find(T == Y_new);
                            end

                        case 'grid'
                            % grid search
                            % add random offset
                            if iteration == 1
                                grid_def = round(1:length(T)/grid_number:length(T));
                                grid_random_offset = randi(length(T) - grid_def(end)+1)-1;
                                grid_random_startlocation = randi(grid_number);
                            end % if isempty(grid_random_offset)
                            X_new = grid_def(mod(length(X_ind) + grid_random_startlocation, grid_number)+1);

                        case 'random'
                            % random sampling
                            X_new = randi(length(T));

                        case 'EI'
                            % Choose the next stimulation parameter with the expected improvement.
                            [X_new,acq] = tms_nextpoint_ei(y, y_var, GP_fit.postmeanX);

                        case 'probEI'
                            % Choose the next stimulation parameter with the expected improvement.
                            [~,acq] = tms_nextpoint_ei(y, y_var, GP_fit.postmeanX);

                            if sum(isnan(acq))
                                X_new = randi(length(T));
                                acq = NaN(size(T));
                                disp('selecting random simulation target due to nonconvergence')
                            else

                                % select the next target randomly from the acq,
                                % interpreting the acq as the probability of
                                % the target being selected

                                Y_new = randsrc(1,1,[T ; (acq/sum(acq))']);
                                X_new = find(T == Y_new);
                            end


                    end

                    %% Save data
                    [~, temp] = max(y);
                    simulation_results_opt_estimate(repetition_current, ...
                        iteration, acquisition_idx, model_idx, subject_idx) = T(temp);
                    [~, temp2] = max(ground_truth.mean);
                    noncirc_error = abs(T(temp) - ground_truth.X(temp2));
                    circ_error = min(noncirc_error, range(T)-noncirc_error);
                    simulation_results_error(repetition_current, ...
                        iteration, acquisition_idx, model_idx, subject_idx) = circ_error;
                    
                    %% ONLINE PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % must be deleted to allow parallel processing

                    % if settings.plotting
                    % 
                    %     % check if plot exists, else build it
                    %     if ~sum(ismember(findall(0,'type','figure'),online_plot))
                    %         online_plot = figure('Name','Measured Data',...
                    %             'NumberTitle','off');
                    %     end
                    % 
                    %     subplot(2,1,1)
                    %     hold off
                    %     % plot data points
                    %     scatter(T(X_ind), t, 'filled', 'k', 'MarkerFaceAlpha', 0.2, 'Markeredgealpha', 0.2)
                    %     xlim([-pi pi])
                    %     hold on
                    % 
                    %     % plot ground truth
                    %     plot(ground_truth.X, ground_truth.mean, '--k')
                    %     Y_gt = [ground_truth.mean+ground_truth.sd, flip(ground_truth.mean-ground_truth.sd)];
                    %     X_area = [ground_truth.X, flip(ground_truth.X)];
                    %     fill(X_area, Y_gt, 'm', 'LineStyle', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.1)
                    % 
                    %     % plot estimate
                    %     plot(T,y, 'Color', [160 20 20]/255, 'LineWidth',2)
                    %     Y = [(y+y_var)', flip(y-y_var)'];
                    %     fill(X_area, Y, 'm', 'LineStyle', 'none', 'FaceColor', [160 20 20]/255, 'FaceAlpha', 0.1)
                    % 
                    %     % plot next data point
                    %     xline(T(X_new), 'LineWidth',2)
                    % 
                    %     % plot acquisition function if converged
                    %     acq_plot = acq - min(acq);
                    %     acq_plot = acq_plot / max(acq_plot);
                    %     try 
                    %         plot(T,(acq_plot - 1 + min(Y_gt)), ':k', 'LineWidth',2);
                    %         no_acq = 0;
                    %     catch 
                    %         no_acq = 1;
                    %     end 
                    % 
                    %     title(['Simulation fitting ' current_model ' with acq ' current_acquisition ': ' num2str(length(t)) ' samples'])
                    %     xlabel('Phase [rad]')
                    %     if no_acq
                    %         ylim([min(Y_gt), max(Y_gt)])
                    %     else
                    %         ylim([min(Y_gt) - abs(min(acq_plot-1)), max(Y_gt)])
                    %     end
                    %     grid on; drawnow
                    % 
                    %     % plot error over iterations
                    %     subplot(2,1,2)
                    %     hold off
                    %     measured_values = simulation_results_error ...
                    %         (:, :, acquisition_idx, model_idx, subject_idx);
                    %     plot(measured_values', 'LineWidth',1, 'Color', [0 0 0 0.2])
                    %     hold on
                    %     plot(mean(measured_values, 'omitnan'), '-o', ...
                    %         'LineWidth',2, 'Color', [160 20 20]/255, 'MarkerFaceColor',[160 20 20]/255)
                    % 
                    %     grid on
                    %     ylabel('error score [rad]')
                    %     xlabel('iteration')
                    %     xlim([1 max(sum(~isnan(measured_values(1,:))),2)])
                    %     ylim([0 pi])
                    % 
                    % end % if settings.plotting
                end % for iteration = 1:settings.n_iterations
            end % for repetition_current = 1:settings.repetitions
            clear i temp Y temp temp2 circ_error noncirc_error m_N S_N beta phi_test t_new X_new

            % Intermediate save
            mkdir('\Results\')
            save(['\Results\Intermediate_Results_' current_datetime], '-v7.3');
        end % for acquisition_idx = 1:length(settings.type_acquisition)
    end % for model_idx = 1:length(settings.type_model)
end % subject_idx = 1:length(settings.subjects) % iterate over subject datasets

%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
close all

disp('start saving data')
save([settings.path_save '\Results\simulation_output'], '-v7.3');
disp('saving completed.')



