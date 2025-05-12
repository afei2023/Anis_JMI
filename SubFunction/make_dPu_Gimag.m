function dPu0=make_dPu_Gimag(P_down,G_imag,velocity, velocity_nmo, velocity_horizontal, ...
        propagation_taper, num_x, dx, num_z,omega, dz, num_freq, min_freq, max_freq, num_freq_bins)

% Initialize kx (spatial frequency in x-direction)
num_x2 = 2^nextpow2(num_x); % Ensure num_x2 is a power of 2 for FFT efficiency
dkx = 2 * pi / (num_x2 * dx); % Frequency resolution in x-direction
kx = [0:dkx:num_x2 / 2 * dkx, -(num_x2 / 2 - 1) * dkx:dkx:-dkx]'; % Spatial frequencies in x-direction

% read out number of sources
nsrc=size(P_down,3);
nrcv=size(P_down,2);

dPu0=zeros(nrcv,nsrc,num_freq); % Define surface pertubation
prev_wavefield=zeros(nrcv,nsrc,num_freq_bins); % define extrabolator
for iz=num_z:-1:1
    % Extract down-going wavefield for current depth level and frequency range
    temp_wavefield = squeeze(P_down(iz + 1, :, :, min_freq:max_freq));

    Gr_at_depth = G_imag(:, iz + 1); % reflectivity for that depth lvl
    GR_at_depth=repmat(Gr_at_depth,[1 nsrc num_freq_bins]); % repeat for all freq and src
    prev_wavefield=prev_wavefield+GR_at_depth.*temp_wavefield; % add previous wavefield ... total wavefiled
    

    % Velocity models for current depth level
    velocity_at_depth = velocity(:, iz).'; % Velocity estimate for this level
    velocity_nmo_at_depth = velocity_nmo(:, iz).'; % NMO velocity for this level
    velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % Horizontal velocity for this level

    prev_wavefield=For_Propagation(prev_wavefield, velocity_at_depth, velocity_nmo_at_depth, velocity_horizontal_at_depth, ...
        propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);
end
dPu0(:,:,min_freq:max_freq)=prev_wavefield; % save pretubation at surface

