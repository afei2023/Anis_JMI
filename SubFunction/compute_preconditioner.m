function preconditioner = compute_preconditioner(P_down, min_freq_index, max_freq_index, ...
                                                max_freq, velocity_estimate, dx, dz)

% Inputs:
%   P_down            : Downward-propagating wavefield 
%   min_freq_index    : Minimum frequency index
%   max_freq_index    : Maximum frequency index
%   max_freq          : Maximum frequency value (Hz)
%   velocity_estimate : Estimated velocity model for filter design
%   dx, dz            : Grid spacing in x and z directions
%
% Outputs:
%   preconditioner         : Preconditioning weights (to multiply with imaging gradient)
%   sensitivity_smoothed   : Smoothed sensitivity map

    % Step 1: Compute total energy (sensitivity map)
    % Sum over frequency and source dimensions
    sensitivity_map = sum(sum(abs(P_down(:,:,:,min_freq_index:max_freq_index)).^2, 4), 3);

    % Step 2: Apply Gaussian smoothing to the sensitivity map
    cutoff_freq = max_freq / 2; % Gaussian filter cutoff frequency
    sensitivity_smoothed = smooth_gradient(sensitivity_map, cutoff_freq, velocity_estimate, dx, dz);

    % Step 3: Compute preconditioner
    epsilon = 0.25 * max(sensitivity_smoothed(:)); % Stabilization term to avoid division by very small numbers
    preconditioner = 1 ./ (sensitivity_smoothed.' + epsilon); % Note: transpose to match gradient dimension

end
