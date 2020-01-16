function stat_energ_plot(grid_size, def_wavenum, lay1_file_name, lay2_file_name)
% Generates a plot of statistical energy versus wavelength, as shown in 
% Figure 5 of Di Qi's paper. Given the disturbance potential voriticity 
% in the barotropic and baroclinic modes, denoted pot_vort_trop and 
% pot_vort_clin, we calculate their corresponding streamfunctions in 
% Fourier space, denoted strmfunc_trop and strmfunc_clin, and compute the 
% arrays containing the barotropic and baroclinic statistical energies 
% (given by the expressions in the caption of Figure 5). We then take the 
% radial average of these matrices, and plot them.

% Frequency-space operators used to obtain the streamfunctions for the 
% disturbance in potential vorticity in the barotropic and baroclinic
% modes.
wavenumbers = [0:grid_size/2 -grid_size/2+1:-1]';
[x_wavenumbers, y_wavenumbers] = meshgrid(wavenumbers, wavenumbers);
freq_deriv_x = 1i*repmat(wavenumbers',[grid_size 1 2]);
freq_deriv_y = 1i*repmat(wavenumbers,[1 grid_size 2]);
freq_laplacian = freq_deriv_x(:,:,1).^2+freq_deriv_y(:,:,1).^2;
inv_freq_trop = 1./freq_laplacian; inv_freq_trop(1,1) = 0;
inv_freq_clin = 1./(freq_laplacian-def_wavenum^2); inv_freq_clin(1,1) = 0;

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

%strmfunc_trop = fftshift(strmfunc_trop);
%strmfunc_clin = fftshift(strmfunc_clin);

% Calculate the radial average of the statisical energy of the barotropic
% and baroclinic modes.
stat_energ_trop = zeros(grid_size/2+1,1);
stat_energ_clin = zeros(grid_size/2+1,1);
for i = 1:grid_size
    for j = 1:grid_size
        wavenum = sqrt(x_wavenumbers(i,j)^2 + y_wavenumbers(i,j)^2);
        if ceil(wavenum) <= grid_size/2
            radius_bin = wavenum - floor(wavenum);
            stat_energ_trop(floor(wavenum)+1) = stat_energ_trop(floor(wavenum)+1) ...
                + (1-radius_bin)*(wavenum^2)*abs(strmfunc_trop(i,j))^2;
            stat_energ_trop(ceil(wavenum)+1) = stat_energ_trop(ceil(wavenum)+1) ...
                + radius_bin*(wavenum^2)*abs(strmfunc_trop(i,j))^2;
            stat_energ_clin(floor(wavenum)+1) = stat_energ_clin(floor(wavenum)+1)...
                + (1-radius_bin)*(wavenum^2 + def_wavenum^2)*abs(strmfunc_clin(i,j))^2;
            stat_energ_clin(ceil(wavenum)+1) = stat_energ_clin(ceil(wavenum)+1)...
                + radius_bin*(wavenum^2 + def_wavenum^2)*abs(strmfunc_clin(i,j))^2;
        end
    end
end

stat_energ_trop = 0.5*stat_energ_trop/(grid_size^4);
stat_energ_clin = 0.5*stat_energ_clin/(grid_size^4);

% Plot the energy matrix for the barotropic part.
%h = pcolor(stat_energ_trop(:,:));
%set(h, 'EdgeColor', 'none');
%disp(stat_energ_trop(1:5,1:5));

% Plot the radial average for the barotropic part.
figure;
loglog(stat_energ_trop,'.-');

% Compute the radial average power spectrum of each statistical energy matrix.
%[Zr_trop,R_trop] = radialavg(stat_energ_trop,grid_size/4,0,0);
%figure;loglog(R_trop,Zr_trop,'.-');

%raPsd2d(stat_energ_trop,0.0078125);
