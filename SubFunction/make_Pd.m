function P_down = make_Pd(source_matrix, P_up, reflectivity, velocity, velocity_nmo, velocity_horizontal, ...
    propagation_taper, num_x, dx, num_z, omega, dz, freq_min_idx, freq_max_idx,num_freq_bins)
    % Inputs:
    % - source_matrix: Source matrix (S2)
    % - P_up: Up-going wavefield
    % - reflectivity: Reflectivity model (imag2)
    % - velocity: Velocity model (vel2)
    % - velocity_nmo: NMO velocity model (vel_nmo)
    % - velocity_horizontal: Horizontal velocity model (vel_h)
    % - propagation_taper: Propagation taper function (Wtaper)
    % - num_x, dx: Grid size and spacing in x-direction
    % - num_z, dz: Grid size and spacing in z-direction
    % - num_time_samples, time_interval: Number of time samples and time interval
    % - freq_min_idx, freq_max_idx: Frequency range indices (ifmin, ifmax)
    %
    % Outputs:
    % - P_down: Down-going incident wavefield (P_down)

    % Initialize kx (spatial frequency in x-direction)
    num_x2 = 2^nextpow2(num_x); % Ensure num_x2 is a power of 2 for FFT efficiency
    dkx = 2 * pi / (num_x2 * dx); % Frequency resolution in x-direction
    kx = [0:dkx:num_x2 / 2 * dkx, -(num_x2 / 2 - 1) * dkx:dkx:-dkx]'; % Spatial frequencies in x-direction

    % Preallocate wavefields
    num_sources = size(source_matrix, 2);
    P_down = zeros(size(P_up)); % Down-going wavefield
    previous_wavefield = zeros(num_x, num_sources, num_freq_bins); % Temporary storage for wavefield

    % Loop through depth levels (from top to bottom)
    for iz = 1:num_z
        % Extract up-going wavefield for current depth level and frequency range
        temp_wavefield = squeeze(P_up(iz, :, :, freq_min_idx:freq_max_idx));

        % Reflectivity for current depth level
        reflectivity_at_depth = reflectivity(:, iz);
        reflectivity_matrix = repmat(reflectivity_at_depth, [1, num_sources, num_freq_bins]);
        % Update the previous wavefield based on current data
        if iz == 1
            previous_wavefield = (1 + reflectivity_matrix) .* previous_wavefield - reflectivity_matrix .* temp_wavefield + source_matrix(:, :, freq_min_idx:freq_max_idx); % Initial scatterer + source
        else
            previous_wavefield = (1 + reflectivity_matrix) .* previous_wavefield - reflectivity_matrix .* temp_wavefield; % Secondary scatterers
        end

        % Velocity models for current depth level
        velocity_at_depth = velocity(:, iz).'; % Velocity for this level
        velocity_nmo_at_depth = velocity_nmo(:, iz).'; % NMO velocity for this level
        velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % Horizontal velocity for this level

        % Propagate the wavefield using the velocity model and the propagation taper
        previous_wavefield = For_Propagation(previous_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
            velocity_horizontal_at_depth, propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);

        % Store the down-going wavefield (P_down)
        P_down(iz + 1, :, :, freq_min_idx:freq_max_idx) = previous_wavefield; % Save updated wavefield
    end
end

