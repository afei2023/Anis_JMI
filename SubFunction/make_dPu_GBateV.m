function dPu0=make_dPu_GBateV(Q_up, G_Beta, velocity, velocity_nmo, velocity_horizontal, ...
        propagation_taper, num_x, dx, num_z,omega, dz, num_freq, min_freq, max_freq,num_freq_bins)

% Initialize kx (spatial frequency in x-direction)
num_x2 = 2^nextpow2(num_x); % Ensure num_x2 is a power of 2 for FFT efficiency
dkx = 2 * pi / (num_x2 * dx); % Frequency resolution in x-direction
kx = [0:dkx:num_x2 / 2 * dkx, -(num_x2 / 2 - 1) * dkx:dkx:-dkx]'; % Spatial frequencies in x-direction

% read out number of sources
nsrc=size(Q_up,3);
nrcv=size(Q_up,2); 

dPu0=zeros(nrcv,nsrc,num_freq);
prev_wavefield=zeros(nrcv,nsrc,num_freq_bins); % define extrabolator
for iz=num_z:-1:1
    qu=squeeze(Q_up(iz+1,:,:,min_freq:max_freq)); % Q- for that depth lvl and freq
    
    velocity_at_depth = velocity(:, iz).'; % Velocity estimate for this level
    velocity_nmo_at_depth = velocity_nmo(:, iz).'; % NMO velocity for this level
    velocity_horizontal_at_depth = velocity_horizontal(:, iz).'; % Horizontal velocity for this level

    % extrapolate  
    prev_wavefield = For_Propagation(prev_wavefield, velocity_at_depth, velocity_nmo_at_depth, ...
       velocity_horizontal_at_depth, propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);

    % apply B scattering
    gbate=G_Beta(:,iz); % beta 
    Gbate=repmat(gbate,[1 nsrc num_freq_bins]); 

    dpu = For_GV_operator(qu, velocity_at_depth, velocity_nmo_at_depth, velocity_horizontal_at_depth, ...
        propagation_taper, omega, kx, num_x, num_x2, dx, dz, num_freq_bins);

    dpu=Gbate.*dpu; % beta*GWQmv
     
    prev_wavefield=prev_wavefield+dpu;  % add to the next step lvl  
end
dPu0(:,:,min_freq:max_freq)=prev_wavefield;