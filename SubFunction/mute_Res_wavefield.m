function Res_wavefield = mute_Res_wavefield(Res_wavefield, min_freq_index, max_freq_index, taper_width, nx)
% mute_wavefield
% Directly mute wavefield outside desired frequency and spatial regions
%
% Inputs:
%   Res_wavefield   - 3D wavefield (x, z, frequency)
%   min_freq_index  - minimum frequency index to preserve
%   max_freq_index  - maximum frequency index to preserve
%   taper_width     - number of muted grid points at each lateral boundary
%   nx              - number of original x grid points (without taper)
%
% Output:
%   Res_wavefield   - muted wavefield

    % Mute frequency components outside the desired band
    Res_wavefield(:, :, 1:min_freq_index-1) = 0;
    Res_wavefield(:, :, max_freq_index+1:end) = 0;

    % Mute spatial regions beyond taper boundaries
    Res_wavefield(1:taper_width, :, :) = 0;
    Res_wavefield(taper_width + nx + 1:end, :, :) = 0;

end
