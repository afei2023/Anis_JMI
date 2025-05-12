function Output = Back_GH_operator(Res_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
    velocity_horizontal_at_depth, propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins)
% Inputs:
%   Res_wavefield                 : Residual wavefield at current step
%   velocity_at_depth             : Velocity profile at current depth
%   velocity_nmo_at_depth         : NMO velocity profile at current depth
%   velocity_horizontal_at_depth  : Horizontal velocity profile at current depth
%   propagation_taper             : Spatial tapering operator
%   omega                         : Frequency vector (angular frequency)
%   kx                            : Wavenumber vector (spatial frequency in x)
%   num_x                         : Number of grid points in x
%   num_x2                        : Next power-of-2 number of x for FFT
%   dx                            : Grid spacing in x-direction
%   dz                            : Grid spacing in z-direction
%   num_freq_bins                 : Number of frequency bins

%% Step 1: Prepare and Preprocessing

% Basic quantities
slowness = 1 ./ velocity_at_depth;
slowness_squared = slowness.^2;
x_grid = (1:num_x) * dx - dx; % Spatial x-coordinates

% Anisotropic parameters
eta = (velocity_horizontal_at_depth.^2 - velocity_nmo_at_depth.^2) ./ (2 * velocity_nmo_at_depth.^2);
delta = 0.5 * ((velocity_nmo_at_depth.^2 ./ velocity_at_depth.^2) - 1);

% Expand into 3D
omega_squared_3D = permute(repmat(omega.^2, [num_x2, 1, num_x]), [1, 3, 2]);
kx_3D = repmat(kx, [1, num_x, num_freq_bins]);
x_grid_3D = repmat(x_grid, [num_x2, 1, num_freq_bins]);
slowness_squared_3D = repmat(slowness_squared, [num_x2, 1, num_freq_bins]);
eta_3D = repmat(eta, [num_x2, 1, num_freq_bins]);
delta_3D = repmat(delta, [num_x2, 1, num_freq_bins]);
velocity_nmo_3D = repmat(velocity_nmo_at_depth, [num_x2, 1, num_freq_bins]);
velocity_horizontal_3D = repmat(velocity_horizontal_at_depth, [num_x2, 1, num_freq_bins]);
velocity_3D = repmat(velocity_at_depth, [num_x2, 1, num_freq_bins]);

%% Step 2: Build Propagation Gradient Operator G

% Dispersion relation terms
A_ps = omega_squared_3D - velocity_nmo_3D.^2 .* (1 + 2 * eta_3D) .* kx_3D.^2;
B_ps = omega_squared_3D - 2 * velocity_nmo_3D.^2 .* eta_3D .* kx_3D.^2;
C_ps = omega_squared_3D .* slowness_squared_3D;

% kz (vertical wavenumber) with anisotropy
kz_3D = conj(sqrt(C_ps .* (A_ps ./ B_ps)));

% Stabilization term
stab_eps = 0.01;
epsilon = stab_eps * omega_squared_3D .* slowness_squared_3D;

% Anisotropic correction term D_ani
D_ani = (omega_squared_3D .* velocity_3D.^2 .* velocity_nmo_3D.^2 .* velocity_horizontal_3D.^2 .* kx_3D.^4 .* conj(kz_3D)) ./ ...
       ((omega_squared_3D .* velocity_3D.^2 - velocity_3D.^2 .* kx_3D.^2 .* (velocity_horizontal_3D.^2 - velocity_nmo_3D.^2)).^2 .* (abs(kz_3D).^2 + epsilon));

% Gradient operator G in frequency domain
Gkx = -1i * 0.5 * dz * (D_ani) .* exp(-1i * kz_3D * dz) .* exp(-1i * kx_3D .* x_grid_3D);

% Remove pseudo S-wave modes
Gkx(A_ps < 0 & B_ps < 0) = 0;

% ifft to spatial domain
Gx = ifft(Gkx);
Gx = Gx(1:num_x, :, :) .* repmat(propagation_taper, [1, 1, num_freq_bins]);

%% Step 3: Apply Forward Operator G to Input Wavefield

Output = zeros(size(Res_wavefield));

parfor ifreq = 1:num_freq_bins
    wavefield_slice = squeeze(Res_wavefield(:, :, ifreq));
    gradient_operator = squeeze(Gx(:, :, ifreq));

    % Apply adjoint operator G
    Output(:, :, ifreq) = gradient_operator * wavefield_slice;
end

end

