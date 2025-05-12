function wave = wavelet(f0, nt, dt)
% Generate Ricker wavelet in frequency domain
% f0: peak frequency (Hz)
% nt: number of time samples
% dt: time sampling interval (s)

nf = nt/2 + 1; % number of positive frequencies
df = 1 / (nt * dt); % frequency interval

% Frequency vector
f = (0:nf-1) * df;

% Ricker wavelet spectrum (positive frequencies)
wave = 2 * f.^2 ./ (sqrt(pi) * f0^3) .* exp(-f.^2 / f0^2);

% Optional: to get full frequency domain [0, ..., f_Nyquist, -f_Nyquist, ..., -df]
% wave_full = [wave, conj(wave(end-1:-1:2))];
end