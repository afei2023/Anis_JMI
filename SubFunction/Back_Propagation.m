function Pout = Back_Propagation(previous_wavefield, velocity_at_depth, velocity_nmo_at_depth, velocity_horizontal_at_depth, ...
    propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins)

% Function to back propagate wavefields
% Inputs:
%   previous_wavefield         : Wavefield at current step
%   velocity_at_depth          : Velocity profile at current depth
%   velocity_nmo_at_depth      : NMO velocity profile at current depth
%   velocity_horizontal_at_depth : Horizontal velocity profile
%   propagation_taper          : Spatial tapering operator
%   omega                      : Frequency vector (angular frequency)
%   kx                         : Wavenumber vector (spatial frequency in x)
%   num_x                      : Number of grid points in x
%   num_x2                     : Next power-of-2 number of x for FFT
%   dx                         : Grid spacing in x-direction
%   dz                         : Grid spacing in z-direction
%   num_freq_bins              : Number of frequency bins

% Calculate slowness
slowness = 1 ./ velocity_at_depth;
slowness_squared = slowness.^2;
x_grid = (1:num_x) * dx - dx; % Set up spatial grid

% Calculate Thomsen’s anisotropic parameter eta
eta = (velocity_horizontal_at_depth.^2 - velocity_nmo_at_depth.^2) ./ (2 * velocity_nmo_at_depth.^2);

% Prepare 3D matrices for calculation
omega_squared_3D = permute(repmat(omega.^2, [num_x2, 1, num_x]), [1, 3, 2]);
kx_3D = repmat(kx, [1, num_x, num_freq_bins]);
x_grid_3D = repmat(x_grid, [num_x2, 1, num_freq_bins]);
slowness_squared_3D = repmat(slowness_squared, [num_x2, 1, num_freq_bins]);
eta_3D = repmat(eta, [num_x2, 1, num_freq_bins]);

% Prepare velocity matrices
velocity_nmo_3D = repmat(velocity_nmo_at_depth, [num_x2, 1, num_freq_bins]);

% Calculate anisotropic propagation terms
A_ps = omega_squared_3D - velocity_nmo_3D.^2 .* (1 + 2 * eta_3D) .* kx_3D.^2;
B_ps = omega_squared_3D - 2 * velocity_nmo_3D.^2 .* eta_3D .* kx_3D.^2;
C_ps = omega_squared_3D .* slowness_squared_3D;
kz_3D = conj(sqrt(C_ps .* (A_ps ./ B_ps)));

% Build propagation operator in kx and kz
Wkx = exp(-1i * kz_3D * dz) .* exp(-1i * kx_3D .* x_grid_3D);

% Transform operator to spatial domain (x-space)
Wx = ifft(Wkx);
Wx = Wx(1:num_x, :, :); % Crop to physical domain

% Apply tapering to suppress boundary artifacts
taper_3D = repmat(propagation_taper, [1, 1, num_freq_bins]);
Wx = taper_3D .* Wx;

% Initialize output
Pout = zeros(size(previous_wavefield));

% Apply propagation operator to each frequency
parfor ifr = 1:num_freq_bins
    wavefield_slice = squeeze(previous_wavefield(:,:,ifr)); % Extract wavefield at frequency ifr
    propagation_operator = squeeze(Wx(:,:,ifr));            % Extract propagation operator at frequency ifr

    % Apply propagation (adjoint operator W')
    Pout(:,:,ifr) = propagation_operator' * wavefield_slice;
end

end

