function S = make_S(wavelet, sxloc, nx, nt, dt)

nsrc = length(sxloc);        % number of sources
nf = nt/2 + 1;               % number of positive frequencies
dom = 2*pi / (nt*dt);        % frequency step in radians

% Initialize source matrix
S = zeros(nx, nsrc, nf);

% Place 1 at each source location for all frequencies
for isrc = 1:nsrc
    S(sxloc(isrc), isrc, :) = 1;
end

% Set DC (0 Hz) and Nyquist components to zero
S(:,:,1) = 0;
S(:,:,nf) = 0;

% Frequency vector (skip DC and Nyquist)
freq_indices = 2:nf-1;
omega = (freq_indices - 1) * dom; % rad/s

% Apply wavelet and frequency scaling
for k = 1:length(freq_indices)
    ifr = freq_indices(k);
    S(:,:,ifr) = wavelet(ifr) * S(:,:,ifr) / sqrt(1i * omega(k));
end

end
