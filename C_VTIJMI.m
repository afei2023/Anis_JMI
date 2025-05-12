%% model migration and inversion 
clc;
clear;
close all;

addpath('Subroutines')
addpath('NSST_denoise')
parpool(8);

%% Simulation Parameters
num_roundtrips = 2;             % Number of scattering round trips
min_freq = 5;                   % Minimum frequency in Hz
max_freq1 = 20;
max_freq2 = 30;
max_freq3 = 40;                 % Maximum frequency in Hz 

dt = 0.004;                      % Time sampling interval (seconds)
dx = 20;                         % Grid spacing in x-direction (meters)
dz = 10;                         % Grid spacing in z-direction (meters)

frequency_band_iterations = [24, 20, 16];  % Iterations for different frequency bands [5~20] [5~30] [5~40]
special_iterations_fraction = 4;  % After 1/4 of iterations in each band for updating Vh

%% Load Models and Data
load Data/FWdata.mat            % load Observed data
load Data/FWwave.mat            % load same ricker wave
load Data/cmap.mat;             % Load colormap
load Data/cgray.mat;            % Load gray colormap

load('Data/vel.mat');           % Load velocity model
load('Data/del.mat');           % Load delta model
load('Data/epsilon.mat');       % Load epsilon model
load('Data/imag.mat');          % Load reflectivity model

%% Model Dimensions and Grid Definitions
num_x = size(vel, 1);           % Number of grid points in x-direction
num_z = size(vel, 2);           % Number of grid points in z-direction
num_freq = size(P_obs, 3);    % Number of frequency
num_sources = size(P_obs, 2); % Number of sources

nt = 2 * num_freq - 2; % number of time samples
freq_interval = 1 / (nt * dt); % Frequency interval
domega = 2 * pi / (nt * dt);   % Angular frequency increment

% Generate spatial and time vectors
x_grid = 0 : dx : (num_x - 1) / 2 * dx;                     % X-location vector
time_grid = 0 : dt : (nt - 1) * dt;                         % Time vector
depth_grid = 0 : dz : (num_z - 1) * dz;                     % Depth vector

%% Ricker Wavelet Source matrix Generation
source_locations = round(linspace(1, num_x, num_sources));  % Define source locations
source_matrix = make_S(ricker_wavelet, source_locations, num_x, nt, dt); % Generate source matrix

%% Extend Models with Tapering
angle_max = 80;                    % Maximum angle for taper
taper_width = 30;                   % Taper width at model boundaries

extended_num_x = num_x + 2 * taper_width;   % Extended model size in x-direction
[extended_vel, extended_delta, extended_epsilon, extended_imag, extended_source_matrix] = ...
    extend_models_with_tapering(vel, delta, epsilon, imag, source_matrix, taper_width, num_x, num_sources, num_freq);
extended_P_obs = zeros(extended_num_x, num_sources, num_freq); % Initialize Data_obs
extended_P_obs(taper_width + 1 : taper_width + num_x, :, :) = P_obs; % Fill in tapered Data_obs

velocity_vertical = extended_vel;                           % ture vertical velocity
velocity_nmo = velocity_vertical .* sqrt(1 + 2 * extended_delta); % ture NMO velocity model !!!Constant
velocity_horizontal = velocity_vertical .* sqrt(1 + 2 * extended_epsilon); % ture horizontal velocity model

% Create angle taper for propagation
propagation_taper = create_angle_taper(angle_max, taper_width, num_x, dx, dz);

%% Set flags for updating Vh %%
total_iterations = sum(frequency_band_iterations);

% Preallocate arrays for frequency and special iteration flags
max_freq_per_iter = zeros(1, total_iterations);  % Store max frequency for each iteration
iter_flag = zeros(1, total_iterations);  % Flags for Vh update

% Get cumulative sum for the end index of each frequency band
end_idx = cumsum(frequency_band_iterations);  % Cumulative sum gives the end index for each band
start_idx = [1, end_idx(1:end-1) + 1];  % Start index for each band is the end index of the previous band + 1

% Assign max frequency for each iteration without using a loop
max_freq_per_iter(1:end_idx(1)) = max_freq1;  % First band
max_freq_per_iter(end_idx(1)+1:end_idx(2)) = max_freq2;  % Second band
max_freq_per_iter(end_idx(2)+1:end_idx(3)) = max_freq3;  % Third band

% Assign flags for each band
iter_flag(1:round(end_idx(1)/special_iterations_fraction)) = 1;  % First band
iter_flag(end_idx(1)+1:round((end_idx(2)-end_idx(1))/special_iterations_fraction) + end_idx(1)) = 1;  % Second band
iter_flag(end_idx(2)+1:round((end_idx(3)-end_idx(2))/special_iterations_fraction) + end_idx(2)) = 1;  % Third band

%% Normalise data %%
% Ensure all frequencies contribute equally in inversion
% Calculate the weighting function
spec = abs(ricker_wavelet(1:num_freq)).^2;
weight = sqrt(spec) ./ (spec + 0.1 * max(spec)); 

% Apply the weight to source and observed data
extended_source_matrix = extended_source_matrix .* reshape(weight, 1, 1, []);
extended_P_obs = extended_P_obs .* reshape(weight, 1, 1, []);

%% Initial Velocity Model (Custom Definition)
velocity_vertical0 = zeros(size(velocity_vertical));
velocity_vertical0(:, 1:15) = 1800;  % First part of velocity model
velocity_vertical0(:, 16:end) = repmat(linspace(1800, 2400, num_z - 15), [extended_num_x, 1]); % Second part

velocity_horizontal0 = zeros(size(velocity_horizontal));  
velocity_horizontal0(:, 1:15) = 1885;     % First part of horizontal velocity model
velocity_horizontal0(:, 16:end) = repmat(linspace(1885, 2550, num_z - 15), [extended_num_x, 1]); % Second part

%%  Prepare loop parameter
% initialise image and velocities
imag_estimate = zeros(size(extended_imag));
velocity_estimate = velocity_vertical0; % inital Vv with taper
velocityh_estimate = velocity_horizontal0;% % initial Vh with taper

% initialise contrast functions
BetaV = zeros(size(velocity_estimate)); % beta 
Betah = zeros(size(velocity_estimate));
P_up = zeros(num_z+1,extended_num_x,num_sources,num_freq); % P-
ObjFun = zeros(1,total_iterations); %

%% define size of figures %%
% define screen resolution
screen_width = 1920;
screen_height = 1080;

% Window size and offset parameters
window_division = 3;

% Define window positions and sizes (x, y, width, height)
window_positions = [
    0, screen_height/2, screen_width/window_division, screen_height/2;                % figure 1
    screen_width/window_division, screen_height/2, screen_width/window_division, screen_height/2;   % figure 2
    2*screen_width/window_division, screen_height/2, screen_width/window_division, screen_height/2; % figure 3
    2*screen_width/3, 0, screen_width/5, screen_height/3;           % figure 4
    0, 0, screen_width/2, screen_height/3;                        % figure 5    
];

% Create all figure windows
for i = 1:size(window_positions, 1)
    figure(i);
    set(gcf, 'position', floor(window_positions(i, :)));
end

drawnow;

for iter=1:total_iterations % loop over number of itrations
    iter
    tic
    % set frequency range and angular frequency
    max_freq=max_freq_per_iter(iter); % max freq to use in this itration
    min_freq_index = round(min_freq / freq_interval) + 1; % Min frequency index--constant
    max_freq_index = round(max_freq / freq_interval) + 1; % Max frequency index
    num_freq_bins = max_freq_index - min_freq_index + 1;  % Number of frequency bins
    omega = domega * (min_freq_index:max_freq_index) - domega;  % Angular frequency vector

    % update vel_nmo
    velocity_nmo = velocity_vertical .* sqrt(1 + 2 * extended_delta);
    %% Imaging
    
    % Update wavefields ... output size nt ... extra freq are padded
    P_down = make_Pd(extended_source_matrix, P_up, imag_estimate, velocity_estimate, velocity_nmo, velocityh_estimate, propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);
    P_up = make_Pu(P_down, imag_estimate, velocity_estimate, velocity_nmo, velocityh_estimate, propagation_taper, extended_num_x, dx, num_z, omega, dz, min_freq_index, max_freq_index,num_freq_bins);
    
    % extract upgoing wavefield at the surface
    surface_wavefield = squeeze(P_up(1, :, :, :));
    % calculate residual
    Res = extended_P_obs - surface_wavefield;
    % mute outside frequencies and edges
    Res_wavefeild = mute_Res_wavefield(Res, min_freq_index, max_freq_index, taper_width, num_x);

    % calculate image gradient
    G_imag = make_Gimag(Res_wavefeild,P_down,velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);

    %  Preconditioning the reflectivity gradient === Guassian smoother
    preconditioner = compute_preconditioner( ...
                                P_down, min_freq_index, max_freq_index, ...
                                max_freq, velocity_estimate, dx, dz);
    G_imag = G_imag .* preconditioner; % Apply Preconditioner to Gradient

    % calculate wavefield perturbation due to G_image
    dPu0 = make_dPu_Gimag(P_down,G_imag,velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, num_freq, min_freq_index, max_freq_index, num_freq_bins);
    taper_indices = [1:taper_width, taper_width+num_x+1:extended_num_x];
    dPu0(taper_indices, :, :) = 0;
        
    % calculate steplength for update image
    dPu0_freq = squeeze(dPu0(:,:,min_freq_index:max_freq_index));
    Res_wavefield_freq = squeeze(Res_wavefeild(:,:,min_freq_index:max_freq_index));
    teller = sum(sum(sum(real(conj(dPu0_freq) .* Res_wavefield_freq))));
    noemer = sum(sum(sum(abs(dPu0_freq).^2)));
    alpha_imag = teller / noemer;

    % update image
    imag_estimate = imag_estimate + alpha_imag * G_imag;
    
    % plot image
    figure(1); clf
    colormap(cgray);
    cmax = max(abs(imag(:)));
    subplot(2,1,1); imagesc(x_grid,(0:num_z)*dz,imag.',[-cmax cmax]); 
      title(['reflectivity, iteration ' num2str(iter)]); xlabel('lateral location [m]'); ylabel('depth [m]');
    cmax = max(abs(imag_estimate(:)));
    subplot(2,1,2); imagesc(x_grid,(0:num_z)*dz,imag_estimate(taper_width+1:taper_width+num_x,:).',[-cmax cmax]); 
      title('Estimated Reflectivity'); xlabel('lateral location [m]'); ylabel('depth [m]');
    drawnow;
    pause(0.1)

    %% Velocity updating
    % Update wavefields
    P_down = make_Pd(extended_source_matrix, P_up, imag_estimate, velocity_estimate, velocity_nmo, velocityh_estimate, propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);
    [P_up, Q_up] = make_Pu(P_down, imag_estimate, velocity_estimate, velocity_nmo, velocityh_estimate, propagation_taper, extended_num_x, dx, num_z, omega, dz, min_freq_index, max_freq_index,num_freq_bins);
  
    % extract upgoing wavefield at the surface
    surface_wavefield = squeeze(P_up(1, :, :, :));
    % calculate residual
    Res = extended_P_obs - surface_wavefield;
    % mute outside frequencies and edges
    Res_wavefeild = mute_Res_wavefield(Res, min_freq_index, max_freq_index, taper_width, num_x);

    % calculate BetaV gradient
    G_BetaV = make_GBetaV(Res_wavefeild,Q_up,velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);

    %  Preconditioning the BetaV(vertical velocity) gradient === Guassian smoother
    preconditioner = compute_preconditioner( ...
                                Q_up(2:num_z+1,:,:,:), min_freq_index, max_freq_index, ...
                                max_freq, velocity_estimate, dx, dz);
    G_BetaV = G_BetaV .* preconditioner; % Apply Preconditioner to Gradient
    G_BetaV = smooth_gradient(G_BetaV, max_freq / 2,velocity_estimate / 2, dx, dz); % Gaussain smoother for G_BetaV
    
    % calculate wavefield perturbation due to slowness update in Pmin0
    dPu0_Vv=make_dPu_GBateV(Q_up, G_BetaV, velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, num_freq, min_freq_index, max_freq_index,num_freq_bins);
    taper_indices = [1:taper_width, taper_width+num_x+1:extended_num_x];
    dPu0_Vv(taper_indices, :, :) = 0;
        
    % calculate steplength for update image
    dPu0_Vv_freq = squeeze(dPu0_Vv(:,:,min_freq_index:max_freq_index));
    Res_wavefield_freq = squeeze(Res_wavefeild(:,:,min_freq_index:max_freq_index));
    teller = sum(sum(sum(real(conj(dPu0_Vv_freq) .* Res_wavefield_freq))));
    noemer = sum(sum(sum(abs(dPu0_Vv_freq).^2)));
    alpha_BetaV = teller / noemer;

    %%%%%%%%%%%%%%%%%%%%%%% Vh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iter_flag(iter) == 0

    % calculate BetaV gradient
    G_BetaH = make_GBetaH(Res_wavefeild,Q_up,velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);

    %  Preconditioning the BetaV(vertical velocity) gradient === Guassian smoother
    preconditioner = compute_preconditioner( ...
                                Q_up(2:num_z+1,:,:,:), min_freq_index, max_freq_index, ...
                                max_freq, velocityh_estimate, dx, dz);
    G_BetaH = G_BetaH .* preconditioner; % Apply Preconditioner to Gradient
    G_BetaH = smooth_gradient(G_BetaH, max_freq / 3, velocity_estimate / 2, dx, dz); % Gaussain smoother for G_BetaV
    
    % calculate wavefield perturbation due to slowness update in Pmin0
    dPu0_Vh=make_dPu_GBateH(Q_up, G_BetaH, velocity_estimate, velocity_nmo, velocityh_estimate, ...
        propagation_taper, extended_num_x, dx, num_z,omega, dz, num_freq, min_freq_index, max_freq_index,num_freq_bins);
    taper_indices = [1:taper_width, taper_width+num_x+1:extended_num_x];
    dPu0_Vh(taper_indices, :, :) = 0;

    % Calculate update step size
    Vv_sub = dPu0_Vv(:, :, min_freq_index:max_freq_index);
    Vh_sub = dPu0_Vh(:, :, min_freq_index:max_freq_index);
    Res_sub = Res_wavefeild(:, :, min_freq_index:max_freq_index);

    Vv_energy = sum(abs(Vv_sub).^2, 'all');
    Vh_energy = sum(abs(Vh_sub).^2, 'all');
    VhVv_cross = sum(real(Vh_sub .* Vv_sub), 'all');

    PdR_vv = sum(real(conj(Vv_sub) .* Res_sub), 'all');
    PdR_vh = sum(real(conj(Vh_sub) .* Res_sub), 'all');

    alpha = [Vv_energy, VhVv_cross; VhVv_cross, Vh_energy]\ [PdR_vv; PdR_vh];
    alpha_betaV = alpha(1);
    alpha_betaH = alpha(2);

    % update Vv and Vh
    velocityh_estimate = velocityh_estimate ./ sqrt(1 + alpha_betaH .* G_BetaH);
    velocity_estimate  = velocity_estimate  ./ sqrt(1 + alpha_betaV .* G_BetaV);

else
   
   Vv_sub = dPu0_Vv(:, :, min_freq_index:max_freq_index);
   Res_sub = Res_wavefeild(:, :, min_freq_index:max_freq_index);
   teller = sum(real(conj(Vv_sub) .* Res_sub), 'all');
   noemer = sum(abs(Vv_sub).^2, 'all');
   alpha_betaV=teller/noemer;
   % update Vv
   velocity_estimate = velocity_estimate ./ sqrt(1 + alpha_betaV .* G_BetaV);

end

    %% Plotting
    figure(2); colormap('default');
    vmin=min(min([velocity_vertical velocity_horizontal]))-300;
    vmax=max(max([velocity_vertical velocity_horizontal]))+300;
    subplot(2,1,1); imagesc(x_grid,(0:num_z)*dz,velocity_vertical(taper_width+1:taper_width+num_x,:).',[vmin vmax]); 
    title(['true vertical velocity, iteration ' num2str(iter)]); xlabel('lateral location [m]'); ylabel('depth [m]');
    subplot(2,1,2); imagesc(x_grid,(0:num_z)*dz,velocity_estimate(taper_width+1:taper_width+num_x,:).',[vmin vmax]);
    title(['inverted vertical velocity, iteration ' num2str(iter)]); xlabel('lateral location [m]'); ylabel('depth [m]');

    figure(3); colormap('default');
    vmin=min(min([velocity_vertical velocity_horizontal]))-300;
    vmax=max(max([velocity_vertical velocity_horizontal]))+300;
    subplot(2,1,1); imagesc(x_grid,(0:num_z)*dz,velocity_horizontal(taper_width+1:taper_width+num_x,:).',[vmin vmax]); 
    title(['true horizontal velocity, iteration ' num2str(iter)]); xlabel('lateral location [m]'); ylabel('depth [m]');
    subplot(2,1,2); imagesc(x_grid,(0:num_z)*dz,velocityh_estimate(taper_width+1:taper_width+num_x,:).',[vmin vmax]);
    title(['inverted horizontal velocity, iteration ' num2str(iter)]); xlabel('lateral location [m]'); ylabel('depth [m]');

    pplot=zeros(num_x,nt);
    pplot(:,min_freq_index:max_freq_index)=squeeze(extended_P_obs(taper_width+1:taper_width+num_x,ceil(num_sources/2),min_freq_index:max_freq_index));
    pplot(:,num_freq+1:nt)=conj(pplot(:,num_freq-1:-1:2));
    pplot=real(ifft(pplot,[],2));
    cmax=max(max(abs(pplot)))/2;
 
    figure(4); colormap(cgray);
    pplot=zeros(num_x,nt);
    pplot(:,min_freq_index:max_freq_index)=squeeze(Res_wavefeild(taper_width+1:taper_width+num_x,ceil(num_sources/2),min_freq_index:max_freq_index));
    pplot(:,num_freq+1:nt)=conj(pplot(:,num_freq-1:-1:2));
    pplot=real(ifft(pplot,[],2));
    imagesc(x_grid,time_grid,pplot.',[-cmax cmax]); xlabel('receiver location [m]'); ylabel('time [s]'); title('residual')
    axis square ; drawnow;

    % Cost function
    ObjFun(iter)=sum(sum(sum(abs(Res_wavefeild).^2)))./sum(sum(sum(abs(extended_P_obs).^2)));
    figure(5);    
    plot(1:total_iterations,log10(ObjFun)); title(['Log10(ObjFun) up to iteration ',num2str(iter)]);
    set(gca,'fontsize',20); xlabel('iteration'); ylabel('Log(error)')   

    pause(0.1);

    toc

end

delete(gcp('nocreate'))
