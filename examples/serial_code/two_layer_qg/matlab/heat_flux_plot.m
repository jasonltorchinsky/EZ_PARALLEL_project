function heat_flux_data = heat_flux_plot

% Simulation and output file parameters.
grid_size = 256;
def_wavenum = 4;
num_files = 331;
output_freq = 20;

% Frequency-space operators used to obtain the streamfunctions for the 
% disturbance in potential vorticity in the barotropic and baroclinic
% modes.
wavenumbers = [0:grid_size/2 -grid_size/2+1:-1]';
freq_deriv_x = 1i*repmat(wavenumbers',[grid_size 1 2]);
freq_deriv_y = 1i*repmat(wavenumbers,[1 grid_size 2]);
freq_laplacian = freq_deriv_x(:,:,1).^2+freq_deriv_y(:,:,1).^2;
inv_freq_trop = 1./freq_laplacian; inv_freq_trop(1,1) = 0;
inv_freq_clin = 1./(freq_laplacian-def_wavenum^2); inv_freq_clin(1,1) = 0;
wavenumbers = [0:grid_size/2-1 0 -grid_size/2+1:-1]';
freq_deriv_x = 1i*repmat(wavenumbers',[grid_size 1]);
%freq_deriv_y = 1i*repmat(wavenumbers,[1 grid_size]);

% Create arrays for storing the time and heat flux data.
time_data = zeros(num_files,1);
heat_flux_data = zeros(num_files,1);

% Loop through all available output files.
for file_num = 0:num_files-1
    lay1_file_name = sprintf('out1_%08d.csv', output_freq*file_num);
    lay2_file_name = sprintf('out2_%08d.csv', output_freq*file_num);
    timestep_data_file_name = sprintf('timestep_info_%08d.csv', output_freq*file_num);
    
    
    % Read the files containing the disturbance in potential vorticity in the
    % barotropic and baroclinic modes.
    pot_vort_lay1 = dlmread(lay1_file_name);
    pot_vort_lay2 = dlmread(lay2_file_name);
    pot_vort_lay1 = pot_vort_lay1(:,1:end-1);
    pot_vort_lay2 = pot_vort_lay2(:,1:end-1);
    
    pot_vort_trop = 0.5 * (pot_vort_lay1 + pot_vort_lay2);
    pot_vort_clin = 0.5 * (pot_vort_lay1 - pot_vort_lay2);
    
    % Calculate the corresponding steamfunctions in Fourier Space.
    pot_vort_trop = fft2(pot_vort_trop);
    pot_vort_clin = fft2(pot_vort_clin);
    
    strmfunc_trop = inv_freq_trop.*pot_vort_trop;
    strmfunc_clin = inv_freq_clin.*pot_vort_clin;
    
    % Calculate the real-space quantities.
    merid_vel = real(ifft2(freq_deriv_x.*strmfunc_trop));
    real_strmfunc_clin = real(ifft2(strmfunc_clin));
    
    % Calculate the area-integrated meridional heat flux.
    heat_flux = ((2 * pi * def_wavenum/grid_size)^2) ...
        *sum(sum(merid_vel.*real_strmfunc_clin));
    
    % Get current time from timestep file.
    [label, data] = readvars(timestep_data_file_name);
    time = data(1);
    
    % Store times heat flux data in their respective arrays.
    time_data(file_num+1) = time;
    heat_flux_data(file_num+1) = heat_flux;
    
    
end


% Plot the time series evolution of heat flux.
plot(time_data,heat_flux_data, 'k-')
hold on
title("Time-Series of Heat Flux (Low-/Mid-Latitude Atmosphere)")
xlabel("Time")
ylabel("Heat Flux")

hold off