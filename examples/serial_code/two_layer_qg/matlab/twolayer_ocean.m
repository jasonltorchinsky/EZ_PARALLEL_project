function twolayer_ocean
% This script solves nondimensional 2-layer QG with equal layers and a
% rigid lid in a doubly-periodic domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set simulation parameters
N = 256;     % Number of points in each direction
dt = 0.0005;   % initial time step size
Nt = 6600;      % Number of time steps
qlim = 150000.;  % if any q > qlim, simulation stops

% Set physical parameters
kd = 4;     % Nondimensional deformation wavenumber
kb = sqrt(2.5);      % Nondimensional beta wavenumber, beta = kb^2
U = 0.2;
r = 0.05;       % Nondimensional Ekman friction coefficient
nu = 5*10^(-15);    % Coefficient of biharmonic vorticity diffusion

outputfreq = 20;

% Set up hyperviscous PV dissipation
k = [0:N/2 -N/2+1:-1]';  % wavenumbers
L = zeros([N N 2]);
for jj=1:N
    for ii=1:N
        kr = sqrt(k(ii)^2+k(jj)^2);
        L(ii,jj,:) = -nu*kr^8;
    end
end
clear kr ii jj


% Put useful stuff into a struct
params = struct('U',U, 'kd',kd, 'kb',kb, 'r',r, 'nu',nu, 'N',N, 'dt',dt);

% Initialize
t = 0;
%qp(:,:,2) = 10*randn(params.N);
%qp(:,:,2) = qp(:,:,2)-mean(mean(qp(:,:,2)));
x0 = 0.0;
y0 = 0.0;
dx = 0.0078125;
dy = 0.0078125;
for kk = 1:N
    for jj = 1:N
        x = x0 + (jj-1)*dx;
        y = y0 + (kk-1)*dy;
        qp(kk,jj,1) = cos(2*x*pi) * sin(2*y*pi);
        qp(kk,jj,2) = cos(2*x*pi) * sin(2*y*pi);
    end
end
q = fft2(qp);

% filename1 corresponds to layer 1, filename2 corresponds to layer 2,
% filename3 corresponds to timestep information.
filename1 = strcat('out1_', sprintf('%08d', 0), '.csv');
filename2 = strcat('out2_', sprintf('%08d', 0), '.csv');
filename3 = strcat('timestep_info_', sprintf('%08d', 0), '.csv');
fileID1 = fopen(filename1, 'w');
fileID2 = fopen(filename2, 'w');
fileID3 = fopen(filename3, 'w');
for kk = 1:N
  for jj = 1:N
    fprintf(fileID1,' %.16f', qp(kk,jj,1));
    fprintf(fileID1,' %s', ', ');
    fprintf(fileID2,' %.16f', qp(kk,jj,2));
    fprintf(fileID2,' %s', ', ');
  end
  fprintf(fileID1, '\n');
  fprintf(fileID2, '\n');
end

fprintf(fileID3, 'Current time: ');
fprintf(fileID3,' %s', ', ');
fprintf(fileID3, ' %.16f', t);
fprintf(fileID3, '\n');
fprintf(fileID3, 'Current time step size: ');
fprintf(fileID3,' %s', ', ');
fprintf(fileID3, ' %.16f', dt);
fprintf(fileID3, '\n');
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);

% Diagnostics
% countDiag = 2; % Compute diagnostics every countDiag steps
% T = zeros(1,Nt/countDiag);
% vb = zeros(1,Nt/countDiag);
% utz = zeros(N,Nt/countDiag);
% ke = zeros(N/2+1,Nt/countDiag);
% ape = zeros(N/2+1,Nt/countDiag);


% adaptive stepping stuff:
tol= 1E-1;
r0 = .8*tol;
% Main loop 
tic;
ii = 0;
while ii < Nt
    ii = ii + 1;
    disp(['Beginning step ' num2str(ii) '.']);
%     if mod(ii,countDiag)==0
%         if any(isnan(q(:))), break, end
%         T(ii/countDiag)=t;
%         [KE,APE] = Spectrum(q,params);
%         ke(:,ii/countDiag) = KE; ape(:,ii/countDiag) = APE;
%         [VB,UTZ] = QG_Diagnostics(q,params);
%         vb(ii/countDiag) = VB; utz(:,ii/countDiag) = UTZ;
%         if mod(ii, countDiag)==0
%             display(['iteration i = ', num2str(ii), '; time step dt = ',num2str(dt), ', ene = ',num2str(sum(KE+APE))]);
%         end
%         
%     end
    M = 1./(1-.25*dt*L);
    % First stage ARK4
    k0 = RHS_Spectral(q,params);
    l0 = L.*q;
    % Second stage
    q1 = M.*(q+.5*dt*k0+.25*dt*l0);
    k1 = RHS_Spectral(q1,params);
    l1 = L.*q1;
    % Third stage
    q2 = M.*(q+dt*(13861*k0/62500+6889*k1/62500+8611*l0/62500-1743*l1/31250));
    k2 = RHS_Spectral(q2,params);
    l2 = L.*q2;
    % Fourth stage
    q3 = M.*(q+dt*(-0.04884659515311858*k0-0.1777206523264010*k1+0.8465672474795196*k2...
    +0.1446368660269822*l0-0.2239319076133447*l1+0.4492950415863626*l2));
    k3 = RHS_Spectral(q3,params);
    l3 = L.*q3;
    % Fifth stage
    q4 = M.*(q+dt*(-0.1554168584249155*k0-0.3567050098221991*k1+1.058725879868443*k2...
    +0.3033959883786719*k3+0.09825878328356477*l0-0.5915442428196704*l1...
    +0.8101210538282996*l2+0.2831644057078060*l3));
    k4 = RHS_Spectral(q4,params);
    l4 = L.*q4;
    % Sixth stage
    q5 = M.*(q+dt*(0.2014243506726763*k0+0.008742057842904184*k1+0.1599399570716811*k2...
    +0.4038290605220775*k3+0.2260645738906608*k4+0.1579162951616714*l0...
    +0.1867589405240008*l2+0.6805652953093346*l3-0.2752405309950067*l4));
    k5 = RHS_Spectral(q5,params);
    l5 = L.*q5;
    
    % Error control
    SCRTCH_pre = 0.003204494398459*(k0+l0) -0.002446251136679*(k2+l2)-0.021480075919587*(k3+l3)...
         +0.043946868068572*(k4+l4) -0.023225035410765*(k5+l5);
    SCRTCH = ifft2(SCRTCH_pre);
     r1 = dt*max(max(max(abs(SCRTCH))));
     if r1>tol
         disp(['r1 = ' num2str(r1) ' is too big... Restarting step ' num2str(ii) '.']);
         dt = .75*dt;
         ii = ii-1;
         continue
     end


            
    % Successful step, proceed to evaluation
    t = t+dt;
    disp(['Time at step ' num2str(ii) ' = ' num2str(t) '.']);
    disp(['dt at step ' num2str(ii) ' = ' num2str(dt) '.']);
    disp(['Error at step ' num2str(ii) ' = ' num2str(r1) '.']);
    qp = real(ifft2(q+dt*(0.1579162951616714*(k0+l0)+0.1867589405240008*(k2+l2)+...
    0.6805652953093346*(k3+l3)-0.2752405309950067*(k4+l4)+(k5+l5)/4)));
    q = fft2(qp);
    % step size adjustment: EPS, PI.3.4 ; divide by 4 for a 4th order
    % method with 3rd order embedded
     dt = ((.75*tol/r1)^(.3/4))*((r0/r1)^(.4/4))*dt;
     disp(['dt for step ' num2str(ii+1) ' = ' num2str(dt) '.']);
     % Print time step data to file.
     %filename3 = strcat('outinfo_', sprintf('%08d', ii), '.txt');
     %fileID3 = fopen(filename3, 'w');
     %fclose(fileID3);
     %fprintf(fileID3,' %
     
     r0=r1;
    if any(abs(qp(:))>qlim)
        disp(['qp = ', num2str(max(abs(qp(:)))),'\n']);
        break
    end
    %Work on writing output array to file with high precision.
    if mod(ii,outputfreq) == 0
        disp(['Outputting step ' num2str(ii) '.']);
        % filename1 corresponds to layer 1, filename2 corresponds to layer 2,
        % filename3 corresponds to timestep information.
        filename1 = strcat('out1_', sprintf('%08d', ii), '.csv');
        filename2 = strcat('out2_', sprintf('%08d', ii), '.csv');
        filename3 = strcat('timestep_info_', sprintf('%08d', ii), '.csv');
        fileID1 = fopen(filename1, 'w');
        fileID2 = fopen(filename2, 'w');
        fileID3 = fopen(filename3, 'w');
        for kk = 1:N
            for jj = 1:N
                fprintf(fileID1,' %.16f', qp(kk,jj,1));
                fprintf(fileID1,' %s', ', ');
                fprintf(fileID2,' %.16f', qp(kk,jj,2));
                fprintf(fileID2,' %s', ', ');
            end
            fprintf(fileID1, '\n');
            fprintf(fileID2, '\n');
        end
        
        fprintf(fileID3, 'Current time: ');
        fprintf(fileID3,' %s', ', ');
        fprintf(fileID3, ' %.16f', t);
        fprintf(fileID3, '\n');
        fprintf(fileID3, 'Current time step size: ');
        fprintf(fileID3,' %s', ', ');
        fprintf(fileID3, ' %.16f', dt);
        fprintf(fileID3, '\n');
        fclose(fileID1);
        fclose(fileID2);
        fclose(fileID3);
    end
end

toc;

%if any(isnan(q(:)))
%    fprintf('NaN\n')
%else
%    save(['oce_pert0','N',num2str(N),'U',num2str(U*10),'kb',num2str(floor(kb)),'kd',num2str(floor(kd)),'_r',num2str(r*10)],...
%                'ii','countDiag','dt','tol','params','T','ke','ape','utz','vb','qp');
%end
end