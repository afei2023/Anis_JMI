function [P_up, Q_up] = make_Pu(P_down, reflectivity, velocity, velocity_nmo, velocity_horizontal, ...
    propagation_taper, num_x, dx, num_z, omega, dz, freq_min_idx, freq_max_idx,num_freq_bins)
    % Inputs:
    % - P_down: Down-going wavefield
    % - reflectivity: Reflectivity model
    % - velocity: Estimated velocity model
    % - velocity_nmo: NMO velocity model
    % - velocity_horizontal: Horizontal velocity model
    % - propagation_taper: Propagation taper function
    % - num_x, dx: Grid size and grid spacing in x-direction
    % - num_z, dz: Grid size and grid spacing in z-direction
    % - num_time_samples, time_interval: Number of time samples and time interval
    % - freq_min_idx, freq_max_idx: Frequency range indices
    %
    % Outputs:
    % - P_up: Up-going incident wavefield
    % - Q_up: Up-going outgoing wavefield

    % Initialize kx (spatial frequency in x-direction)
    num_x2 = 2^nextpow2(num_x); % Ensure num_x2 is a power of 2 for FFT efficiency
    dkx = 2 * pi / (num_x2 * dx); % Frequency resolution in x-direction
    kx = [0:dkx:num_x2 / 2 * dkx, -(num_x2 / 2 - 1) * dkx:dkx:-dkx]'; % Spatial frequencies in x-direction

    % Preallocate wavefields
    num_sources = size(P_down, 3);
    P_up = zeros(size(P_down)); % Up-going wavefield
    Q_up = zeros(size(P_down)); % Additional wavefield to store previous wavefield data
    prev_wavefield = zeros(num_x, num_sources, num_freq_bins); % Temporary storage for previous wavefield

    % Loop through depth levels (from bottom to top)
    for iz = num_z:-1:1
        % Extract down-going wavefield for current depth level and frequency range
        temp_wavefield = squeeze(P_down(iz + 1, :, :, freq_min_idx:freq_max_idx));

        % Reflectivity for current depth level
        reflectivity_at_depth = reflectivity(:, iz + 1);
        reflectivity_matrix = repmat(reflectivity_at_depth, [1, num_sources, num_freq_bins]); % Replicate reflectivity for each source and frequency

        % Update the previous wavefield based on current data
        prev_wavefield = (1 - reflectivity_matrix) .* prev_wavefield + reflectivity_matrix .* temp_wavefield; % Add previous wavefield to current wavefield
        Q_up(iz + 1, :, :, freq_min_idx:freq_max_idx) = prev_wavefield; % Store updated wavefield

        % Velocity models for current depth level
        velocity_at_depth = velocity(:, iz).'; % Velocity estimate for this level
        velocity_nmo_at_depth = velocity_nmo(:, iz).'; % NMO velocity for this level
        velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % Horizontal velocity for this level

        % Propagate the wavefield using the velocity model and the propagation taper
        prev_wavefield = For_Propagation(prev_wavefield, velocity_at_depth, velocity_nmo_at_depth, velocity_horizontal_at_depth, ...
            propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);

        % Store the up-going wavefield (P_up)
        P_up(iz, :, :, freq_min_idx:freq_max_idx) = prev_wavefield; % Save updated wavefield
    end
end

