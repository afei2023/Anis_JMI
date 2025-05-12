function grad_filtered = smooth_gradient(grad, fmax, vel, dx, dz)
% SMOOTH_GRADIENT Smooth gradient using Gaussian filter (space domain).
%
% This function smooths the input gradient by applying separable 1D Gaussian
% filters along the x (horizontal) and z (vertical) directions in space domain.
% It is designed for efficient and accurate smoothing in seismic inversion.
%
% Inputs:
%   grad : [nx, nz] 2D matrix, the gradient to be smoothed
%   fmax : scalar, maximum frequency for smoothing (Hz)
%   vel  : [nx, nz] 2D matrix, velocity model (m/s)
%   dx   : scalar, spatial sampling interval in x direction (m)
%   dz   : scalar, spatial sampling interval in z direction (m)
%
% Output:
%   grad_filtered : [nx, nz] 2D matrix, the smoothed gradient
%

    [nx, nz] = size(grad);

    % Find minimum velocity to define maximum wavenumber
    vmin = min(vel(:));

    % Define standard deviations for the Gaussian kernel
    sigma_z = vmin / (2 * fmax * sqrt(2*pi));  % in meters
    sigma_x = vmin / (2 * fmax * sqrt(2*pi));  % in meters

    % Define half filter lengths, covering ~4 sigma
    half_len_z = max(1, round(4 * sigma_z / dz));
    half_len_x = max(1, round(4 * sigma_x / dx));

    % Generate 1D Gaussian filters
    z = (-half_len_z:half_len_z) * dz;
    filterz = exp(-(z.^2) / (2 * sigma_z^2));
    filterz = filterz / sum(filterz);  % normalize to preserve energy

    x = (-half_len_x:half_len_x) * dx;
    filterx = exp(-(x.^2) / (2 * sigma_x^2));
    filterx = filterx / sum(filterx);  % normalize

    % Extend the gradient at the boundaries by replication
    grad_ext = padarray(grad, [half_len_x, half_len_z], 'replicate', 'both');

    % Apply Gaussian smoothing
    grad_tmp = conv2(grad_ext, filterz(:)', 'same');  % smooth along z (columns)
    grad_tmp = conv2(grad_tmp, filterx(:), 'same');   % smooth along x (rows)

    % Crop back to original size
    grad_filtered = grad_tmp(half_len_x+1 : half_len_x+nx, half_len_z+1 : half_len_z+nz);

end




