function init = BO_settings()
% function init = BO_settings()
% input dialogue for settings of BO simulation
% 
% INPUTS
    %
% OUTPUTS
    % init      struct with fields containing settings
        % initial_samples       



prompt = {...
    'Number of initial random samples',...
    'Number of iterations'...
    'Model type(s)'...
    'Acquisition function(s)'...
    'Plotting true/false'...
    'Path for loading data'...
    'Path for saving data'...
    'Kernel type'...
    'Subjects (space-separated)'...
    'Repetitions of each model'...
    };
definput = {...
    '10',...
    '1000'...
    'GP, bayesian regression'...
    'KG, random'...
    'true'...
    '\BO_simulation\Summary_data_allsubjects_complete.mat'...
    '\Bayesian optimization simulation\Results'...
    'periodicMK'...
    '1 7 19 22 32'...
    '1000'...
    };
dlgtitle = 'Settings for bayesian optimization script';
dims = [1 50];

answer = inputdlg(prompt,dlgtitle,dims,definput);
 
% extract answers into response struct
init.n_init_samples     = str2double(answer{1});
init.n_iterations       = str2double(answer{2});
init.plotting           = eval(answer{5});
init.path_load          = answer{6};
init.path_save          = answer{7};
init.type_kernel        = answer{8};
init.subjects           = str2num(answer{9});
init.n_repetitions      = str2double(answer{10});

init.type_model = {};
if contains(answer{3}, 'GP','IgnoreCase',true)
    init.type_model = [init.type_model, {'GP'}]; end
if contains(answer{3}, 'simple regression','IgnoreCase',true)
    init.type_model = [init.type_model, {'simple regression'}]; end
if contains(answer{3}, 'bayesian regression','IgnoreCase',true)
    init.type_model = [init.type_model, {'bayesian regression'}]; end

init.type_acquisition   = {};
if contains(answer{4}, 'KG','IgnoreCase',true)
    init.type_acquisition = [init.type_acquisition, {'KG'}]; end
if contains(answer{4}, 'grid','IgnoreCase',true)
    init.type_acquisition = [init.type_acquisition, {'grid'}]; end
if contains(answer{4}, 'random','IgnoreCase',true)
    init.type_acquisition = [init.type_acquisition, {'random'}]; end
if contains(answer{4}, 'EI','IgnoreCase',true)
    init.type_acquisition = [init.type_acquisition, {'EI'}]; end


init.formatSpec = '%03.0f';




