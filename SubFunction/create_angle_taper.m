function Wtaper = create_angle_taper(angle, nxtap, nx, dx, dz)
    % Function to create an angle taper for the propagation operators
    % This function generates a tapering function based on the maximum angle for the propagation operator.

    % Calculate the number of grid points to taper in the x-direction based on the angle
    nxflat = floor(dz * tan(angle * pi / 180) / dx); % Taper width in x-direction based on angle
    ntap = floor((nx - 2 * nxflat) / 2); % Taper width for the interior

    % Total number of grid points in the x-direction, including taper
    nx2 = nx + 2 * nxtap;

    % Create cosine-based taper for the boundaries
    tap = 0.5 + 0.5 * cos(-pi * (1:ntap) / ntap); % Cosine taper function for the boundary

    % Initial tapering vector
    wtaper = [ones(1, nxflat), tap];

    % Initialize Wtaper as an identity matrix
    Wtaper = diag(ones(1, nx2));

    % Apply the taper to the diagonals of Wtaper
    for idiag = 1:length(wtaper)
        Wtaper = Wtaper + wtaper(idiag) * diag(ones(1, nx2 - idiag), idiag); % Upper diagonal
        Wtaper = Wtaper + wtaper(idiag) * diag(ones(1, nx2 - idiag), -idiag); % Lower diagonal
    end

    % Create second taper for the edges of the model
    tap2 = 0.5 + 0.5 * cos(-pi * (1:nxtap) / nxtap); % Cosine taper for the model edges
    wtaper2 = [tap2(end:-1:1), ones(1, nx), tap2]; % Combine edge taper with middle region

    % Final Wtaper is the product of edge and interior tapers
    Wtaper = diag(wtaper2) * Wtaper * diag(wtaper2); % Apply edge tapering

end
