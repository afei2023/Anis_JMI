function G_imag  = make_Gimag(Res_wavefeild,P_down,velocity, velocity_nmo, velocity_horizontal, ...
    propagation_taper, num_x, dx, num_z,omega, dz, min_freq, max_freq,num_freq_bins)
% make_gradientA_V
% Compute the gradient (imaginary part) for velocity inversion using residual back-propagation.
%
% Inputs:
%   Res_wavefield         - residual wavefield (observed - modeled)
%   P_down                - forward propagated downgoing wavefield
%   velocity              - vertical velocity model
%   velocity_nmo          - NMO velocity model
%   velocity_horizontal   - horizontal velocity model
%   propagation_taper     - taper applied during wave propagation
%   extended_num_x        - number of grid points in x (after model extension)
%   dx                    - spatial sampling interval in x (meters)
%   num_z                 - number of depth levels
%   omega                 - angular frequencies (vector)
%   dz                    - depth sampling interval (meters)
%   min_freq, max_freq    - frequency range to use
%   num_freq_bins         - number of frequency bins
%
% Output:
%   G_imag                - computed gradient (imaginary part)

    % Initialize kx (spatial frequency in x-direction)
    num_x2 = 2^nextpow2(num_x); % Ensure num_x2 is a power of 2 for FFT efficiency
    dkx = 2 * pi / (num_x2 * dx); % Frequency resolution in x-direction
    kx = [0:dkx:num_x2 / 2 * dkx, -(num_x2 / 2 - 1) * dkx:dkx:-dkx]'; % Spatial frequencies in x-direction

   G_imag = zeros(num_x,num_z+1); % Initialize gradiant matrix
   previous_wavefield = Res_wavefeild(:,:, min_freq:max_freq); % Temporary storage for previous wavefield

for iz=1:num_z
    pdown=squeeze(P_down(iz+1,:,:,min_freq:max_freq)); % P+(z+1) and required freq
    
    % Extract velocity models at current depth
    velocity_at_depth = velocity(:, iz).'; % Velocity for this level
    velocity_nmo_at_depth = velocity_nmo(:, iz).'; % NMO velocity for this level
    velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % Horizontal velocity for this level

    previous_wavefield = Back_Propagation(previous_wavefield, ...
            velocity_at_depth, velocity_nmo_at_depth, velocity_horizontal_at_depth, ...
            propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);

    % Cross-correlation imaging
    G_imag(:, iz+1) = sum(sum(4 * real(previous_wavefield .* conj(pdown)), 2), 3);
end