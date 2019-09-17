%% Analysis_idRDK.
% access correct directory, set as root, and list subject names / numbers.
root = 'C:\Users\cooperd\Documents\idRDK\Results'; cd(root)
observers = dir(fullfile(root, '*Conv.mat*'));
Blue_values = [146 123 105 107 125 109 122 127 113 114 126 112 138 143 111 126 111 133 106 102]';

%% create file hierarchy w/ accessibility.
for obs_idx = 1:length(observers) % index all 20 observers.
    load(fullfile(root, observers(obs_idx).name))

%% index ONLY run_1 -- BIAS of assessing ambiguous sphere.
    for run_idx = 1
        load(fullfile(root, strrep(observers(obs_idx).name,'Conv.mat',['run_' num2str(run_idx) '.mat']))) % strrep = string replace
        clear frequency_inner frequency_outer correct congruent

        for trial_idx = 1:length(Results.PDir) % index all 120 trial runs.
            if length(Results.PDir{trial_idx}) > 1 % removes any double inputs.
                clean_PDir(trial_idx) = NaN;
            else
                clean_PDir(trial_idx) = Results.PDir{trial_idx};
            end
        end
        % assess probability of direction switch; looks for differences between inputs NOT EQUAL to zero.
        Probability_Switch_single(obs_idx,1) = length(find(diff(clean_PDir) ~= 0))/Results.trial_n;
        % assess directional bias of single sphere tests.
        Bias_single(obs_idx,1) = (length(find(clean_PDir == +1)))/Results.trial_n;
    end

%% index ONLY run_2 -- CORRECTNESS of assessing disambiguated sphere.
    for run_idx = 2
        load(fullfile(root, strrep(observers(obs_idx).name,'Conv.mat',['run_' num2str(run_idx) '.mat'])))
        clear frequency_inner frequency_outer correct congruent

        for trial_idx = 1:length(Results.PDir) % index all 120 trial runs.
            if length(Results.PDir{trial_idx}) > 1
                clean_PDir(trial_idx) = NaN;
            else
                clean_PDir(trial_idx) = Results.PDir{trial_idx};
            end
        end
        % assess accuracy of subjects; omit participants < 0.8.
        Correct_outer(obs_idx,1) = (length(find((clean_PDir - Results.random_vector) == 0)))/Results.trial_n;
    end

%% index run_3-6 -- BIAS of assessing congruency.
    for run_idx = 3:6
        load(fullfile(root, strrep(observers(obs_idx).name,'Conv.mat',['run_' num2str(run_idx) '.mat'])))
        clear frequency_inner frequency_outer correct congruent

        for trial_idx = 1:length(Results.PDir)
            if length(Results.PDir{trial_idx}) > 1
                clean_PDir(trial_idx) = NaN;
            else
                clean_PDir(trial_idx) = Results.PDir{trial_idx};
            end
        end
        % assess congruency.
        Probability_Congruency(obs_idx, run_idx-2) = length(find(clean_PDir == 2))/Results.trial_n;
        % assess probability of direction switch.
        Probability_Switch_double(obs_idx, run_idx-2) = length(find(diff(clean_PDir .* Results.random_vector) ~= 0))/Results.trial_n;
        Bias_double(obs_idx,run_idx-2) = (length(find((clean_PDir .* Results.random_vector) == +2)))/Results.trial_n;
    end

%% index ALL runs -- assessing blue value effects.
    for run_idx = 1:6
        load(fullfile(root, strrep(observers(obs_idx).name,'Conv.mat',['run_' num2str(run_idx) '.mat'])))
        clear frequency_inner frequency_outer correct congruent

        for blue_idx = 1:length(Blue_values)
            Results.BlueValue = Blue_values(obs_idx);
        end
        % assess SOMETHING???
        Blue_vs_PDir(obs_idx) = Results.BlueValue/Results.trial_n; % this is just the trial from building the blue script (not using PDir)...
    end
end

%% Stats + Correlation coefficients.
% probability of phase switch of inner sphere with / without context.
[H,P,CI,STAT] = ttest2((mean(Probability_Switch_single, 2) - 0.5), (mean(Probability_Switch_double,2) - 0.5));

% probability of congruency.
[H,P,CI,STAT] = ttest(mean(Probability_Congruency, 2) - 0.5);

% directional bias of inner and outer sphere.
[H,P,CI,STAT] = ttest(mean(Bias_single, 2) - 0.5);
[H,P,CI,STAT] = ttest(mean(Bias_double, 2) - 0.5);

% CORRELATION phase switch with / without context.
[r, p] = corrcoef(mean(Probability_Switch_single,2), mean(Probability_Switch_double,2));

% CORRELATION blue intensity vs. phase switch with / without context.
[r, p] = corrcoef(Blue_values, mean(Probability_Switch_single,2));
[r, p] = corrcoef(Blue_values, mean(Probability_Switch_double,2));

% CORRELATION blue intensity vs. bias single / bias double / congruency
[r, p] = corrcoef(Blue_values, mean(Bias_single,2));
[r, p] = corrcoef(Blue_values, mean(Bias_double, 2));
[r, p] = corrcoef(Blue_values, mean(Probability_Congruency,2));

%% Final Figures.
% lsline not working...
figure(1)
subplot(2,1,1) % probability of switch of inner sphere with / without context.
boxplot([mean(Probability_Switch_single, 2) mean(Probability_Switch_double,2)], 'Notch', 'on', 'Labels', {'No Context', 'Context'}),...
    ylim([0 1]), ylabel('Probability of Directional Switch')
subplot(2,1,2) % probability of congruency.
boxplot(mean(Probability_Congruency, 2), 'Notch', 'on', 'Labels', {'Congruency'}), ylim([0 1.1]), ylabel('Proportion Congruent Percepts'),...
    hold on, plot([0 2], [0.5 0.5], 'r'), hold off

figure(2)
subplot(1,3,1) % correctness of outer sphere.
boxplot(Correct_outer, 'Notch', 'on', 'Labels', {'Disambiguated'}), ylim([0.80 1.01]), ylabel('Percent Correct')
subplot(1,3,2) % range of blue values.
boxplot(Blue_values, 'Notch', 'on', 'Labels', {'HFP'}), ylim([100 150]), ylabel('Relative Blue Intensity')
subplot(1,3,3) % directional bias of inner and outer sphere.
boxplot([Bias_single mean(Bias_double,2)], 'Notch', 'on', 'Labels', {'Single', 'Double'}), ylim([0 1.1]), ylabel('Bias')

figure(3) % CORRELATION probability switch with / without context
scatter(mean(Probability_Switch_single,2), mean(Probability_Switch_double,2), 'filled'),...
    xlabel('Prob. Phase Switch Single'), ylabel('Prob. Phase Switch Double')

figure(4) % CORRELATION blue value vs. probability switch with / without context
scatter(Blue_values, mean(Probability_Switch_single, 2)), xlabel('Relative Blue Intensity'), ylabel('Prob. Phase Switch'), hold on,...
scatter(Blue_values, mean(Probability_Switch_double, 2), 'r'), legend('Single', 'Double'), hold off

figure(5)
subplot(1,3,1) % CORRELATION blue value vs. bias_single
scatter(Blue_values, mean(Bias_single, 2)), ylim([0 1.1]), ylabel('Bias Single')
subplot(1,3,2) % CORRELATION blue value vs. bias_double
scatter(Blue_values, mean(Bias_double, 2)), ylim([0 1.1]), xlabel('Relative Blue Intensity'), ylabel('Bias Double')
subplot(1,3,3) % CORRELATION blue value vs. probability congruency
scatter(Blue_values, mean(Probability_Congruency, 2), 'r'), ylim([0 1.1]), ylabel('Congruency')
