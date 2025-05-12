%% Acoustic Full Wavefield Simulation
clc;
clear;
close all;

addpath('Subroutines')

% Initialize parallel pool for multi-core processing
parpool(8);

%% Simulation Parameters
num_roundtrips = 2; % Number of scattering round trips
min_freq = 5; % Minimum frequency in Hz
max_freq = 40; % Maximum frequency in Hz
peak_freq = 20; % Peak frequency in Hz

nt = 512; % Number of time samples
dt = 0.004; % Time sampling interval in seconds
dx = 20; % Grid spacing in x-direction (meters)
dz = 10; % Grid spacing in z-direction (meters)

%% Load Models
load Data/cmap.mat; % Preset colormap
load Data/cgray.mat; % Preset colormap

load('Data/vel.mat'); 
load('Data/del.mat'); 
load('Data/epsilon.mat');
load('Data/imag.mat'); 

%% Frequency and Time Calculation
num_frequencies = nt / 2 + 1; % Number of frequency bins
freq_interval = 1 / (nt * dt); % Frequency interval

num_x = size(vel, 1); % Model size in x-direction
num_z = size(vel, 2); % Model size in z-direction
num_sources = round(num_x / 1); % Number of sources (adjustable)

% Frequency bounds
min_freq_index = round(min_freq / freq_interval) + 1; % Min frequency index
max_freq_index = round(max_freq / freq_interval) + 1; % Max frequency index
num_freq_bins = max_freq_index - min_freq_index + 1;

% Generate spatial and time vectors
x_grid = -(num_x - 1) / 2 * dx : dx : (num_x - 1) / 2 * dx; % X-location from negative to positive
time_grid = 0 : dt : (nt - 1) * dt; % Time vector
depth_grid = 0 : dz : (num_z - 1) * dz; % Depth vector

%% Ricker Wavelet
ricker_wavelet = wavelet(peak_freq, nt, dt); % Generate Ricker wavelet in frequency domain
ricker_wavelet = nt * ricker_wavelet ./ (2 * sum(ricker_wavelet)); % Normalize amplitude

% Source matrix generation
source_locations = round(linspace(1, num_x, num_sources)); % Define source locations
source_matrix = make_S(ricker_wavelet, source_locations, num_x, nt, dt); % Generate source matrix

%% Extend Models with Tapering
angle_max = 80; % Maximum angle for W taper
taper_width = 30; % Tapering width at model boundaries
extended_num_x = num_x + 2 * taper_width;
[extended_vel, extended_delta, extended_epsilon, extended_imag, extended_source_matrix] = extend_models_with_tapering(vel, delta, epsilon, imag, source_matrix, taper_width, num_x, num_sources, num_frequencies);
propagation_taper = create_angle_taper(angle_max, taper_width, num_x, dx, dz);

%% Model Parameters
velocity_estimate = extended_vel;
velocity_nmo = velocity_estimate .* sqrt(1 + 2 * extended_delta); % NMO velocity model
velocity_horizontal = velocity_estimate .* sqrt(1 + 2 * extended_epsilon); % Horizontal velocity model

% Calculate angular frequency
domega = 2 * pi / (nt * dt);
omega = domega * (min_freq_index:max_freq_index) - domega;

%% Preallocate Wavefields
P_up = zeros(num_z + 1, extended_num_x, num_sources, num_frequencies); 
forward_wavefield = zeros(num_x, nt); % Initialize forward wavefield
prev_forward_wavefield = forward_wavefield; % Temporary storage wavefield for figure

% Simulation Loop for Each Scattering Order
for scattering_order = 1:num_roundtrips
    disp(['Scattering order: ' num2str(scattering_order)]);

    % Update wavefields for current scattering order
    P_down = make_Pd(extended_source_matrix, P_up, extended_imag, velocity_estimate, velocity_nmo, velocity_horizontal, propagation_taper, extended_num_x, dx, num_z,omega, dz, min_freq_index, max_freq_index,num_freq_bins);
    P_up = make_Pu(P_down, extended_imag, velocity_estimate, velocity_nmo, velocity_horizontal, propagation_taper, extended_num_x, dx, num_z, omega, dz, min_freq_index, max_freq_index,num_freq_bins);

    % Surface wavefield for P-
    surface_wavefield = squeeze(P_up(1, taper_width + 1 : taper_width + num_x, :, :)); 

    % Convert to time domain
    forward_wavefield = convert_to_time_domain(surface_wavefield, num_x, nt, num_frequencies, num_sources);
    
    % Plot the wavefield evolution
    figure(1); colormap(cgray);
    cmax = max(max(abs(forward_wavefield))) / 2;
    subplot(1, 3, 1); imagesc(x_grid, time_grid, prev_forward_wavefield.', [-cmax cmax]); 
    title(['Pup_0^- (' num2str(scattering_order - 1) ')']);
    subplot(1, 3, 2); imagesc(x_grid, time_grid, forward_wavefield.', [-cmax cmax]); 
    title(['Pup_0^- (' num2str(scattering_order) ')']);
    if scattering_order == 1
       subplot(1, 3, 3);
       imagesc(x_grid, time_grid, forward_wavefield.' - prev_forward_wavefield.', [-cmax cmax]);
       title('Primary');
    else
       subplot(1, 3, 3);
       imagesc(x_grid, time_grid, forward_wavefield.' - prev_forward_wavefield.', [-cmax cmax]);
       title('Internal multiples');
    end
    pause(0.1)

    prev_forward_wavefield = forward_wavefield;
end

% Save the results
P_obs = surface_wavefield;
save('Data/FWdata.mat', 'P_obs');
save('Data/FWwave.mat', 'ricker_wavelet');

% Clean up parallel pool
delete(gcp('nocreate'));
