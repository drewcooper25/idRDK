%% Settings_DdRDK (structure, rotation & stimulus duration).

%% Frames.
fixation_interval = 0.8; Results.fixation_interval = fixation_interval;

trial_n = 120; Results.trial_n = trial_n; % trial_n = 120.
rot_per_trial = 0.125; Results.rot_per_trial = rot_per_trial;
est_ifi = 1/60; desired_time_for_stimulus = 3; % seconds
frames_per_rot = desired_time_for_stimulus / (rot_per_trial * est_ifi);
Results.Stimulus.frames_per_rot = frames_per_rot;
stereo_d = floor(frames_per_rot/100); Results.stereo_d = stereo_d;
shift_per_frame = 2*pi/frames_per_rot;

est_exp_duration = (est_ifi * frames_per_rot * rot_per_trial * trial_n +...
    fixation_interval * (trial_n + 1)) / 60;

%% Setup Sphere.
num_dots = 5000;
Results.Stimulus.num_dots = num_dots;
step_size_dots = 2*pi/num_dots;
r = 1; x_c = 0; y_c = 0; z_c = 0; % radius and axis centers.

theta = 0:step_size_dots:2*pi-step_size_dots; % polar angle; elevation (z->).

mu = 0.5; v = 0.03; % mean and variance
sigma = sqrt(v); % standard deviation
X = sigma .* randn(1,length(theta)) + mu;
out_of_bound = [find(X<0) find(X>1)];
X(out_of_bound) = rand(1,length(out_of_bound));
R_random =  X; % redistribution of pole point clusters.
%R_random = rand(1, length(theta)); % dots clustered @ poles.

phi = (pi .* R_random) - (0.5 * pi); % azimuthal angle (x->y).

% initial dot configuration of sphere.
x = x_c + r .* cos(phi) .* cos(theta); % (x axis).
y = y_c + r .* cos(phi) .* sin(theta); % (y axis).
z = z_c + r .* sin(phi); % (z axis).

fraction_of_surface = 0.25; % King's Apple slice cutouts (0.3).
Results.Stimulus.fraction_of_surface = fraction_of_surface;
cutout = find(abs(x) > fraction_of_surface & abs(y) > fraction_of_surface);
x(cutout) = []; y(cutout) = []; z(cutout) = [];
theta(cutout) = []; phi(cutout) = [];

inner_size_factor = 300; Results.Stimulus.inner_size_factor = inner_size_factor;
outer_size_factor = 450; Results.Stimulus.outer_size_factor = outer_size_factor;

%% Overlap + Offset.
n_overlaps = 8; Results.n_overlaps = n_overlaps;
lag_between_overlaps = frames_per_rot / n_overlaps;
overlaps = 1:lag_between_overlaps:frames_per_rot*rot_per_trial; Results.overlaps = overlaps;

minimal_response_time = 0.2; Results.minimal_response_time = minimal_response_time;
minimal_response_frames = round(minimal_response_time/est_ifi);
% estimated ifi instead of true???

offset_frames = round(frames_per_rot/n_overlaps/2); Results.offset_frames = offset_frames;

%% Loop Prep.
% first row of matrix r_xxx containing xxx values
r_theta(1,:) = theta; r_phi(1,:) = phi;
r_x(1,:) =  x; r_y(1,:) =  y; r_z(1,:) =  z;

for idx = 2:frames_per_rot % beginning @ second frame and going until 1400
    r_theta(idx,:) = r_theta(idx - 1,:) + shift_per_frame;
    r_theta(idx,r_theta(idx,:) >= (2 * pi)) = r_theta(idx,r_theta(idx,:) >= (2 * pi)) - (2 * pi);
    r_phi(idx,:) = r_phi(idx - 1,:);

    r_x(idx,:) = x_c + r.*cos(r_phi(idx,:)).*cos(r_theta(idx,:));
    r_y(idx,:) = y_c + r.*cos(r_phi(idx,:)).*sin(r_theta(idx,:));
    r_z(idx,:) = z_c + r.*sin(r_phi(idx,:));
end

%%