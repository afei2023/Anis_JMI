function Output = For_Propagation(previous_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
    velocity_horizontal_at_depth, propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins)
    % Function to propagate wavefields considering anisotropic media
    % Inputs:
    %   previous_wavefield    : Initial wavefield
    %   velocity_at_depth      : Velocity profile at each depth level
    %   velocity_nmo_at_depth  : NMO velocity at each depth level
    %   velocity_horizontal_at_depth : Horizontal velocity at each depth level
    %   propagation_taper      : Tapering operator
    %   omega                  : Frequency vector
    %   kx                      : Wavenumber in x direction
    %   num_x                  : Number of grid points in x direction
    %   num_x2                 : Optimized grid size in x direction (next power of 2)
    %   dx                     : Grid spacing in x direction
    %   dz                     : Grid spacing in z direction
    %   freq_min_idx, freq_max_idx : Frequency range to evaluate
    slowness = 1 ./ velocity_at_depth; % Slowness (inverse of velocity)
    slowness_squared = slowness.^2;
    x_grid = (1:num_x) * dx - dx; % x spacing

    % Calculate eta (anisotropic parameter) in terms of velocity
    eta = (velocity_horizontal_at_depth.^2 - velocity_nmo_at_depth.^2) ./ (2 * velocity_nmo_at_depth.^2);
    
    % Prepare x-grid and frequency matrix
    omega_squared = permute(repmat(omega.^2, [num_x2, 1, num_x]), [1, 3, 2]); % omega squared (3D)
    kx_3D = repmat(kx, [1, num_x, num_freq_bins]); % kx (repeated over all frequencies)
    x_grid_3D = repmat(x_grid, [num_x2, 1, num_freq_bins]); % x positions (3D)
    slowness_squared_3D = repmat(slowness_squared, [num_x2, 1, num_freq_bins]); % slowness squared (3D)
    eta_3D = repmat(eta, [num_x2, 1, num_freq_bins]); % eta (3D)

    % Anisotropic propagation (remove pseudo S-wave)
    velocity_nmo = velocity_nmo_at_depth; % Use NMO velocity for this step
    velocity_nmo_3D = repmat(velocity_nmo, [num_x2, 1, num_freq_bins]); % NMO velocity (3D)
    
    % Calculate kz using anisotropic formula (Page 35 of doctoral thesis)
    A_ps = omega_squared - velocity_nmo_3D.^2 .* (1 + 2 * eta_3D) .* kx_3D.^2;
    B_ps = omega_squared - 2 * velocity_nmo_3D.^2 .* eta_3D .* kx_3D.^2;
    C_ps = omega_squared .* slowness_squared_3D;
    kz_3D = conj(sqrt(C_ps .* (A_ps ./ B_ps))); % Anisotropic kz propagation

    % Build the propagation operator (Wkz)
    Wkx = exp(-1i * kz_3D * dz) .* exp(-1i * kx_3D .* x_grid_3D); % Propagation operator in kx and kz
    Wx=ifft(Wkx); % Apply IFFT

    % Apply taper to the operator
    taper_3D = repmat(propagation_taper, [1, 1, num_freq_bins]); % Tapering operator
    Wx = taper_3D .* Wx(1:num_x, :, :); %  Obtain the final operator (Wx)

    % Initialize the output
    Output = zeros(size(previous_wavefield));
    
    % Propagate wavefield using parallel computation
    parfor ifr = 1:num_freq_bins
        wavefield_freq = squeeze(previous_wavefield(:,:,ifr)); % Extract the wavefield at frequency ifr
        propagation_operator = squeeze(Wx(:,:,ifr));   % Extract the propagation operator at frequency ifr

        % Apply the propagation operator to the wavefield
        Output(:,:,ifr) = propagation_operator * wavefield_freq;
    end
end

