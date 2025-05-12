function forward_wavefield = convert_to_time_domain(surface_wavefield, num_x, nt, num_frequencies, num_sources)
    % Initialize the temporary wavefield
    temp_wavefield = zeros(num_x, nt);
    
    % Assign frequencies from the surface wavefield
    temp_wavefield(:, 1:num_frequencies) = squeeze(surface_wavefield(:, round((num_sources + 1) / 2), :));
    
    % Conjugate symmetry for time domain
    temp_wavefield(:, num_frequencies + 1:nt) = conj(temp_wavefield(:, num_frequencies - 1:-1:2));
    
    % Convert to time domain
    forward_wavefield = real(ifft(temp_wavefield, [], 2));
end