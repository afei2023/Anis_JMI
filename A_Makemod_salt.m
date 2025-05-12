clear; clc;

% Define grid parameters
dx = 20; % Grid spacing in x-direction
dz = 10; % Grid spacing in z-direction
sizex = 3000; % Model width in meters
sizez = 1000; % Model depth in meters

% Define background properties
background_vel = 2000; % m/s
background_eps = 0.075; 
background_del = 0.1;

% Define layer properties as arrays
layers = struct();
layers(1).vel = 1800;  layers(1).eps = 0.05;   layers(1).del = 0.075;
layers(1).x = [0, 3000, 3000, 0, 0];
layers(1).z = [0, 0, 150, 150, 0];

layers(2).vel = 2500;  layers(2).eps = 0.15;   layers(2).del = 0.2;
layers(2).x = [500,700,900,1100,1300,1500,1700,1900,2100,2300,2500,2600 ];
layers(2).z = [420,330,300,290,285,290,300,290,293,296,350,420];

layers(3).vel = 2500;  layers(3).eps = 0.15;   layers(3).del = 0.2;
layers(3).x = [500,700,900,1100,1300,1500,1700,1900,2100,2300,2500,2600 ];
layers(3).z = [420,480,495,475,455,455,455,475,475,465,445,420];

layers(4).vel = 2300;  layers(4).eps = 0.1;    layers(4).del = 0.15;
layers(4).x = [0, 3000, 3000, 0, 0];
layers(4).z = [750, 750, 1000, 1000, 750];

% Generate grid
x_grid = 0:dx:sizex;
z_grid = 0:dz:sizez;
[Z, X] = meshgrid(z_grid, x_grid);  % Swapped Z and X to match the size

% Initialize velocity and epsilon models
vel = background_vel * ones(size(X));
epsilon = background_eps * ones(size(X));
delta = background_del * ones(size(X));

% Populate the velocity, epsilon, and delta models based on layers
for i = 1:length(layers)
    mask = inpolygon(X, Z, layers(i).x, layers(i).z);
    vel(mask) = layers(i).vel;
    epsilon(mask) = layers(i).eps;
    delta(mask) = layers(i).del;
end

% Define background velocity model
vel0 = background_vel * ones(size(vel));
vel0(:,1:15) = 1800;
vel0(:,16:end) = repmat(linspace(1800, 2400, size(vel0, 2) - 15), [size(vel0, 1), 1]);

% Define initial horizontal velocity model
velh0 = 1885 * ones(size(vel));
velh0(:,16:end) = repmat(linspace(1885, 2550, size(velh0, 2) - 15), [size(velh0, 1), 1]);

% Create imaginary part of the velocity model

nx=size(vel,1); % number of x locations
nz=size(vel,2);
imag=zeros(nx,nz+1);
for ix=1:nx
    for iz=2:nz
          imag(ix,iz)=4*(vel(ix,iz)-vel(ix,iz-1))/(vel(ix,iz)+vel(ix,iz-1));
    end
end
%imag = 4 * (vel(:,2:end) - vel(:,1:end-1)) ./ (vel(:,2:end) + vel(:,1:end-1));

% Visualization of models
fig_titles = {'Vertical Velocity Model (m/s)', 'Initial Vertical Velocity Model (m/s)', ...
              'Initial Horizontal Velocity Model (m/s)', 'Horizontal Velocity Model (m/s)', ...
              'Eps Model (g/cm³)', 'Del Model (m/s)', 'Reflectivity Model (m/s)'};
models = {vel, vel0, velh0, velh0, epsilon, delta, imag};

for i = 1:length(models)
    figure;
    imagesc(x_grid, z_grid, models{i}');
    set(gca, 'YDir', 'reverse');
    colorbar;
    xlabel('Horizontal Distance (m)');
    ylabel('Depth (m)');
    title(fig_titles{i});
    colormap('default');
    axis equal tight;
end

% Save data
save('Data/vel.mat', 'vel');
save('Data/vel0.mat', 'vel0');
save('Data/imag.mat', 'imag');
save('Data/del.mat', 'delta');
save('Data/epsilon.mat', 'epsilon');

