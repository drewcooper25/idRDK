%% Copyright (C) 2019 Veith Weilnhammer, CCM.
%% V2.0 Build by Drew Cooper, MedNeuro MS.c.
clear all
close all

%% Intermittent_dRDK Experiment
root_dir = "/Users/user/Documents/CHARITE/Masters Thesis/IdRDK";
BlueValue = '123'; Results.BlueValue = BlueValue;
ObserverName = 'test_002';
SettingsName = 'IdRDK'; % Drew's King's Apple.
which_run = 0;

%%
%% Run 0 (Test + Tryout):
if which_run == 0
    session = ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 1; transition_probability = 0;

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    % analysis (not stored)
    [frequency_inner frequency_outer congruent correct] = get_conventional_data(1, present_inner, present_outer, Results);
end

% figure(),
% bar(length(find(abs(correct(:,2)) == 1))/length(correct)), ylabel('Fraction Correct'), ylim([0 1.1]), xlim([0 2]), hold on
% plot([0 2], [0.5 0.5], 'r')

%%
%% Run 1, Ambiguous:
if which_run == 1
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 0; transition_probability = 0;

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);

    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%
%% Run 2, Disambiguated:
if which_run == 2
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 0; present_outer = 1; transition_probability = 0; % issues with 'correct' input readout, still registers hits as zeros...

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);
    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%
%% Run 3, Combo etc.:
if which_run == 3
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 1; transition_probability = mean(mean(frequency_inner(:,(1))));

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);
    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%
%% Run 4:
if which_run == 4
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 1; transition_probability = mean(mean(frequency_inner(:,(1))));

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);
    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%
%% Run 5:
if which_run == 5
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 1; transition_probability = mean(mean(frequency_inner(:,(1))));

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);
    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%
%% Run 6:
if which_run == 6
    session= ['run_' num2str(which_run)]; run_idx = which_run;

    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end

    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end

    present_inner = 1; present_outer = 1; transition_probability = mean(mean(frequency_inner(:,(1))));

    clear Results
    [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName, BlueValue, session, present_inner, present_outer, transition_probability);

    %% analysis
    [frequency_inner(:,run_idx) frequency_outer(:,run_idx) congruent(:,run_idx) correct(:,run_idx)] = get_conventional_data(run_idx, present_inner, present_outer, Results);
    save([root_dir 'Results/Results_' ObserverName '_Conv.mat'], 'frequency_inner', 'frequency_outer', 'congruent', 'correct')
end

%%