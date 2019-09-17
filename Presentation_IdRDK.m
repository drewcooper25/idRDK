%% Presentation_IdRDK (computer & monitor controls).

function [Results] = Presentation_IdRDK(root_dir, SettingsName, ObserverName,...
    BlueValue, session, present_inner, present_outer, transition_probability)

%% Get Observer Info.
Results.ObserverName = ObserverName;
Results.session=session;

%% Run Settings File.
Settingsfile = [root_dir 'Settings/' 'Settings_' SettingsName '.m'];
run(Settingsfile);

%% Run Experimental Settings.
Results.present_inner = present_inner; Results.present_outer = present_outer;
Results.transition_probability = transition_probability;
percepts = [-1 1];

for trial = 1:trial_n
    template.discrete_steps{trial} = zeros(1,length(overlaps));

    if isnan(transition_probability)
        template.discrete_steps{trial} = repmat(percepts(Randi(2)),1,length(overlaps));

    else
        template.discrete_steps{trial}(1) = percepts(Randi(2));

        for idx = 2 : length(overlaps)
            if randp([transition_probability 1-transition_probability],1,1) == 1;
                template.discrete_steps{trial}(idx) = (-1) * template.discrete_steps{trial}(idx-1);
            else
                template.discrete_steps{trial}(idx) = template.discrete_steps{trial}(idx-1);
            end
        end
    end

    template.discrete{trial} = overlaps(find([1 diff(template.discrete_steps{trial})] ~= 0));
    template.PDir{trial} = template.discrete_steps{trial}(find([1 diff(template.discrete_steps{trial})] ~= 0));
end

Results.template = template;

%% Screen and Keyboard Setup.
AssertOpenGL;

KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
left = KbName('LeftArrow'); right = KbName('RightArrow');
up = KbName('UpArrow'); down = KbName('DownArrow');
RestrictKeysForKbCheck([escapeKey left right up down]);

PsychDefaultSetup(2); % unintentional color intensity boost from Screen('ColorRange',...)?
Screen('Preference', 'SkipSyncTests', 2); % 0 = sync tests, 2 = no tests.
screenid = max(Screen('Screens')); doublebuffer = 1;
[window, windowRect] = PsychImaging('OpenWindow', screenid, 0,...
    [], 32, doublebuffer + 1, [], 128); Results.Monitor.windowRect = windowRect;
[xCenter, yCenter] = RectCenter(windowRect);

mon_width = 36.5; v_dist = 60; % monitor width + viewing distance (cm).
Results.Monitor.mon_width = mon_width; Results.Monitor.v_dist = v_dist;

ifi = Screen('GetFlipInterval', window); Results.Monitor.ifi = ifi;
numSecs = 1; waitframes = 1;
numFrames = round(numSecs / ifi);
fps = Screen('FrameRate', window);
if fps == 0
    fps = 1/ifi;
end

Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Colors + ramp_colors.
black = BlackIndex(screenid); white = WhiteIndex(screenid);
red = [100 0 0]; green = [0 100 0]; blue = [0 0 str2double(BlueValue)]; 
% instead of 'str2num'.

n_frames_ramp = 1;

ramp_red = [(0:red(1) / n_frames_ramp:red(1))' zeros(n_frames_ramp + 1, 2);...
    repmat(red, frames_per_rot * rot_per_trial - 2 * (n_frames_ramp + 1), 1);...
    (red(1):-red(1) / n_frames_ramp:0)' zeros(n_frames_ramp + 1, 2)];

ramp_green = [zeros(n_frames_ramp + 1, 1) (0:green(2) / n_frames_ramp:green(2))' zeros(n_frames_ramp + 1, 1);...
    repmat(green,frames_per_rot * rot_per_trial - 2 * (n_frames_ramp + 1), 1); ...
    zeros(n_frames_ramp + 1, 1) (green(2):-green(2) / n_frames_ramp:0)' zeros(n_frames_ramp + 1, 1)];

ramp_blue = [zeros(n_frames_ramp + 1, 2) (0:blue(3) / n_frames_ramp:blue(3))';...
    repmat(blue, frames_per_rot * rot_per_trial - 2 * (n_frames_ramp + 1),1); ...
    zeros(n_frames_ramp + 1, 2) (blue(3):-blue(3) / n_frames_ramp:0)'];

%% Presentation Flip.
HideCursor;
Priority(MaxPriority(window));
vbl = Screen('Flip', window); % initial screen flip.

%% Dot Dimensions.
ppd = pi * (windowRect(3) - windowRect(1)) / atan(mon_width / v_dist/2) / 360;
Results.Monitor.ppd = ppd; % ppd = pixels per degree.
dot_w = 0.1; Results.dot_w = dot_w;
fix_r = 0.10; Results.fix_r = fix_r;
dot_size_Pix = dot_w * ppd;
dot_type = 2;

%% Setup Fixation Cross.
fix_coord = [xCenter - fix_r * ppd yCenter + fix_r * ppd];
FixCross = [xCenter - 1, yCenter - fix_r * ppd, xCenter + 1,...
    yCenter + fix_r * ppd; xCenter - fix_r * ppd, yCenter - 1,...
    xCenter + fix_r * ppd, yCenter + 1];

%% Randomize Dismabiguous Direction.
random_vector = [ones(1,trial_n/2) -ones(1,trial_n/2)]; % right and left randomization...
random_vector = random_vector(randperm(length(random_vector))); % continued.
Results.random_vector = random_vector;

%% Fixation + Drawing Loop.
Results.SessionStartTime = GetSecs;

for trial = 1:trial_n
    n_switch = 0;

    Results.PDir{trial} = []; Results.PosSwitch{trial} = []; Results.SwitchTime{trial} = [];

    Screen('FillRect', window, uint8(white), FixCross');

    Screen('Flip', window);
    WaitSecs(fixation_interval);
    Results.TrialStartTime{trial} = GetSecs;

    idx = 0;
    Results.offset_frames = offset_frames;
    d1_idx = idx + offset_frames;

    for frames = 1:frames_per_rot*rot_per_trial
        Screen('FillRect', window, ramp_red(frames,:) + ramp_blue(frames,:), FixCross');

        idx = idx + 1; % indexes for angle rotation of inner sphere plot.

        d1_idx = d1_idx + random_vector(trial); % direction changes in stimulus.
        d2_idx = d1_idx + stereo_d;

        if idx > frames_per_rot
            idx = idx - frames_per_rot;
        end

        if d1_idx > frames_per_rot
            d1_idx = d1_idx - frames_per_rot;
        end

        if d2_idx > frames_per_rot
            d2_idx = d2_idx - frames_per_rot;
        end

        if idx <= 0
            idx = idx + frames_per_rot;
        end

        if d1_idx <= 0
            d1_idx = d1_idx + frames_per_rot;
        end

        if d2_idx <= 0
            d2_idx = d2_idx + frames_per_rot;
        end

        if present_outer
            Screen('DrawDots', window, [r_x(d1_idx,:); r_z(d1_idx,:)] .* outer_size_factor,...
                dot_size_Pix/2, ramp_blue(frames,:), [xCenter yCenter], dot_type)
            Screen('DrawDots', window, [r_x(d2_idx,:); r_z(d2_idx,:)] .* outer_size_factor,...
                dot_size_Pix/2, ramp_red(frames,:), [xCenter yCenter], dot_type)
        end

        if present_inner
            Screen('DrawDots', window, [r_x(idx,:); r_z(idx,:)] .* inner_size_factor,...
                dot_size_Pix, ramp_red(frames,:) + ramp_blue(frames,:), [xCenter yCenter], dot_type)
        end

        Screen('DrawingFinished', window); % end of drawing, beginning of data collection.

        [keyIsDown, seconds, keyCode] = KbCheck;
        if keyIsDown
            if length(find(keyCode == 1)) == 1
                if keyCode(escapeKey)
                    break;

                elseif keyCode(left) || keyCode(right) || keyCode(up) || keyCode(down)
                    n_switch = n_switch + 1;
                    Results.PDir{trial}(n_switch) = find(keyCode == 1); % +1: clockwise, -1: counterclockwise, 0: unclear.
                    Results.PosSwitch{trial}(n_switch) = frames;
                    Results.SwitchTime{trial}(n_switch) = seconds;
                end
            end
        end

        if doublebuffer == 1
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
    end

    Results.TrialEndTime{trial} = seconds;

    if isempty(Results.PDir{trial}) % if no response was given.
        Results.PDir{trial} = NaN; % percepts(Randi(2));
    end

    RepDir = diff(Results.PDir{trial});
    Results.PDir{trial}(find(RepDir == 0) + 1) = [];
    Results.PosSwitch{trial}(find(RepDir == 0) + 1) = [];
    Results.SwitchTime{trial}(find(RepDir == 0) + 1) = [];
    % clear response entries.

    Results.PDir{trial}(Results.PDir{trial} == left) = -1; %left
    Results.PDir{trial}(Results.PDir{trial} == right) = +1; %right
    Results.PDir{trial}(Results.PDir{trial} == up) = +2; %up
    Results.PDir{trial}(Results.PDir{trial} == down) = -2; %down
    % rename directions.

    Results.discrete{trial}(1) = 1;
    for idx = 2:length(Results.PosSwitch{trial})
        Results.discrete{trial}(idx) = max(overlaps(overlaps...
            <= Results.PosSwitch{trial}(idx) - minimal_response_frames));
    end

    Results.RT{trial} = (Results.PosSwitch{trial} - Results.discrete{trial}) * ifi;
    Results.discrete_steps{trial} = zeros(1,length(overlaps));

    for idx = 1:length(Results.discrete_steps{trial})
        Results.discrete_steps{trial}(idx) = Results.PDir{trial}...
            (max(find((Results.discrete{trial} <= overlaps(idx)))));
    end
    % round to overlaps.

    [keyIsDown, seconds, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(escapeKey)
            break;
        end
    end
end

%% Final Fixation.
Screen('FillRect', window, uint8(white), FixCross');
Screen('Flip', window);
WaitSecs(fixation_interval);
Results.SessionEndTime = GetSecs;

%% Save Results.
save([root_dir 'Results/Results_' Results.ObserverName '_' num2str(session) '.mat'], 'Results')

%% Back Home.
Priority(0);
ShowCursor
sca;

end

%%