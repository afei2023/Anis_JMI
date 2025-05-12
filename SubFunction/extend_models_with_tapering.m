function [extended_vel, extended_delta, extended_epsilon, extended_imag, extended_source_matrix] = extend_models_with_tapering(vel, delta, epsilon, imag, source_matrix, taper_width, num_x, num_sources, num_frequencies)

    % Extend the model in both x and z directions by tapering the boundary regions
    
    extended_num_x = num_x + 2 * taper_width; % New model size with taper
    extended_vel = [repmat(vel(1,:), taper_width, 1); vel; repmat(vel(end,:), taper_width, 1)]; % Tapered velocity model
    extended_delta = [repmat(delta(1,:), taper_width, 1); delta; repmat(delta(end,:), taper_width, 1)]; % Tapered delta model
    extended_epsilon = [repmat(epsilon(1,:), taper_width, 1); epsilon; repmat(epsilon(end,:), taper_width, 1)]; % Tapered epsilon model
    extended_imag = [repmat(imag(1,:), taper_width, 1); imag; repmat(imag(end,:), taper_width, 1)]; % Tapered reflectivity model

    % Create tapered source matrix
    extended_source_matrix = zeros(extended_num_x, num_sources, num_frequencies); % Initialize tapered source matrix
    extended_source_matrix(taper_width + 1 : taper_width + num_x, :, :) = source_matrix; % Fill in tapered source matrix

end
