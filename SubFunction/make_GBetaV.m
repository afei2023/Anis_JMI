function G_BetaV = make_GBetaV(Res, Q_up, velocity, velocity_nmo, velocity_horizontal, ...
    propagation_taper, num_x, dx, num_z, omega, dz, min_freq, max_freq, num_freq_bins)
% MAKE_GRADIENTB_V: Compute model gradient BetaV using residuals and adjoint fields
%
% Inputs:
%   Res                       : Residual wavefield at the surface
%   Q_up                      : Up-going wavefield from forward modeling
%   velocity                  : Background velocity model
%   velocity_nmo              : NMO velocity model
%   velocity_horizontal       : Horizontal velocity model
%   propagation_taper         : Spatial tapering operator
%   num_x                     : Number of grid points in x
%   dx                        : Spatial grid spacing in x
%   num_z                     : Number of depth levels
%   omega                     : Frequency vector (angular frequency)
%   dz                        : Grid spacing in z
%   min_freq, max_freq         : Frequency indices range (selected frequencies)
%   num_freq_bins             : Number of frequency bins considered

% Set up kx (spatial frequency vector)
num_x2 = 2^nextpow2(num_x); % Next power of 2 for efficient FFT
dkx = 2 * pi / (num_x2 * dx); % Frequency resolution
kx = [0:dkx:num_x2/2 * dkx, -(num_x2/2 - 1) * dkx:dkx:-dkx]'; % Full kx spectrum

% Initialize output
G_BetaV = zeros(num_x, num_z); % Gradient model (BetaV)
Res_wavefield = Res(:, :, min_freq:max_freq); % Select frequencies for residual

% Loop over depth levels
for iz = 1:num_z

    Q_slice = squeeze(Q_up(iz+1, :, :, min_freq:max_freq)); % (x, 1, freq)

    % Extract velocity profiles at current depth
    velocity_at_depth = velocity(:, iz).'; % (1, x)
    velocity_nmo_at_depth = velocity_nmo(:, iz).'; % (1, x)
    velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % (1, x)

    %  Extrapolate residual wavefield downward by one depth step using gradient operator G
    dq_extrap = Back_GV_operator(Res_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
        velocity_horizontal_at_depth, propagation_taper, omega, kx, ...
        num_x, num_x2, dx, dz, num_freq_bins);

    % Cross-correlation to compute gradient
    G_BetaV(:, iz) = sum(sum(4 * real(dq_extrap .* conj(Q_slice)), 2), 3);  % G*E-*Q-=dBate-

    % Update residual wavefield by propagating one step downward
    Res_wavefield = Back_Propagation(Res_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
        velocity_horizontal_at_depth, propagation_taper, omega, kx, ...
        num_x, num_x2, dx, dz, num_freq_bins);

end

end
