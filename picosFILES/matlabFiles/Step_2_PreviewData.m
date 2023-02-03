% Extract data from PICOS++ output:

clear all;
close all;
clc;

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;

% Read MAIN file:
% =========================================================================

%File name:
fileName = '../outputFiles/HDF5/main.h5';

% Get data:
m = HDF2Struct_v2(fileName);

% Number of ranks for fields:
ranksFields = double(m.numMPIsFields);

% Number of ranks for particles:
ranksParticles = double(m.numMPIsParticles);

% Total number of mesh points:
Nx = m.mesh.Nx_IN_SIM;

% Mesh points per MPI:
NxPerMPI = double(m.mesh.Nx_PER_MPI);

% x axis:
x_m  = m.mesh.x_m;

% Number of ions species:
numIonSpecies = m.ions.numberOfParticleSpecies;

% Ion species parameters:
for ii = 1:numIonSpecies
    fieldName = ['species_',num2str(ii)];
    ionParameters{ii} = m.ions.(fieldName);
end

% Get number of time outputs:
% =========================================================================
fileName  = ['../outputFiles/HDF5/FIELDS_FILE_0.h5'];
fields    = HDF2Struct_v2(fileName);
timeSteps = fieldnames(fields);
Nt = length(timeSteps);
timeMax = double(fields.(['t',num2str(Nt-1)]).time);


% Determine the outputted variables:
% =========================================================================
% fields:
vars.fields = fieldnames(fields.t0.fields);

% IONS:
fileName  = ['../outputFiles/HDF5/PARTICLES_FILE_0.h5'];
IONS      = HDF2Struct_v2(fileName);
vars.IONS = fieldnames(IONS.t0.ions.species_1);

vars.output = [vars.fields;vars.IONS];

% Display to terminal:
% =========================================================================
disp('********************************************************************')
disp(['Ranks for fields                 : ', num2str(ranksFields)]);
disp(['Ranks for particles              : ', num2str(ranksParticles)]);
disp(['Number of ION species            : ', num2str(numIonSpecies)]);
disp(['Total simulation time            : ', num2str(timeMax*1e3),' [ms]']);
disp(['Number of time steps outputted   : ', num2str(Nt)]);
disp('********************************************************************')
disp('')
%% Read FIELD files:
% =========================================================================
disp('')
disp('********************************************************************')
disp('Extracting field data...')

% Initialize mesh-defined electromagnetic field variables:
% =======================================================
if find(strcmpi('Bx_m',vars.fields) == 1)
    Bx_m   = zeros(Nx,Nt);
elseif find(strcmpi('dBx_m',vars.fields) == 1)
    dBx_m  = zeros(Nx,Nt);
elseif find(strcmpi('ddBx_m',vars.fields) == 1)
    ddBx_m = zeros(Nx,Nt);
elseif find(strcmpi('ddBx_m',vars.fields) == 1)
    Ex_m   = zeros(Nx,Nt);
end

% Get the data:
% =============
for rr = 1:ranksFields
    
    disp(['Fields rank #', num2str(rr)]);
        
    % File name:
    fileName = ['../outputFiles/HDF5/FIELDS_FILE_',num2str(rr-1),'.h5'];
    
    % Get data:
    fields = HDF2Struct_v2(fileName);
    
    % Data range:
    rng_m = 1 + NxPerMPI*(rr-1) : 1 : NxPerMPI*rr;
    
    % Loop over all time:
    for tt = 1:Nt
        % Obtain time stamp:
        kk = str2double(timeSteps{tt}(2:end)) + 1;
        
        % Mesh-defined fields:
        if sum(strcmpi('Bx_m',vars.output))
            Bx_m(rng_m,kk)   = fields.(timeSteps{tt}).fields.Bx_m.x;
        end
        if sum(strcmpi('dBx_m',vars.output))
            dBx_m(rng_m,kk)  = fields.(timeSteps{tt}).fields.dBx_m.x;
        end
        if sum(strcmpi('ddBx_m',vars.output))
            ddBx_m(rng_m,kk) = fields.(timeSteps{tt}).fields.ddBx_m.x;
        end
        if sum(strcmpi('Ex_m',vars.output))
            Ex_m(rng_m,kk) = fields.(timeSteps{tt}).fields.Ex_m.x;
        end
        
    end
end

disp('Field data extraction completed!')
disp('********************************************************************')
disp('')

%% Read PARTICLE files:
% =========================================================================

disp('')
disp('********************************************************************')
disp('Extracting particle data...')

% Initialize particle data:
% ==========================
for ss = 1:numIonSpecies
    % Number of particles per rank:
    N_CP_MPI = double(ionParameters{1}.N_CP_MPI);

    % Total number of particles:
    N_CP = N_CP_MPI*ranksParticles;
        
    % Initialize data:
    t_p        = zeros(Nt ,1 );
    
    % Particle states:
    if sum(strcmpi('x_p',vars.output))
        x_p{ss}    = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('V_p',vars.output))
        vpar_p{ss} = zeros(N_CP,Nt);
        vper_p{ss} = zeros(N_CP,Nt);  
    end
    
    if sum(strcmpi('a_p',vars.output))
        a_p{ss}    = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('mu_p',vars.output))
        mu_p{ss}   = zeros(N_CP,Nt); 
    end        
        
    % Particle-defined moments:
    if sum(strcmpi('n_p',vars.output))
        n_p{ss}    = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('nv_p',vars.output))
        nv_p{ss}   = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('Tpar_p',vars.output))
        Tpar_p{ss} = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('Tper_p',vars.output))
        Tper_p{ss} = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('Te_p',vars.output))
        Te_p{ss} = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('u_p',vars.output))
        ux_p{ss}   = zeros(N_CP,Nt);
    end
    
    % Mesh-defined moments:
    if sum(strcmpi('n_m',vars.output))
        n_m{ss}    = zeros(Nx,Nt);
    end
    if sum(strcmpi('nv_m',vars.output))
        nv_m{ss}   = zeros(Nx,Nt);
    end
    
    if sum(strcmpi('Tpar_m',vars.output))
        Tpar_m{ss} = zeros(Nx,Nt);
    end
    
    if sum(strcmpi('Tper_m',vars.output))
        Tper_m{ss} = zeros(Nx,Nt);
    end
    
    if sum(strcmpi('Te_m',vars.output))
        Te_m{ss} = zeros(Nx,Nt);
    end
    
    if sum(strcmpi('u_m',vars.output))
        ux_m{ss}   = zeros(Nx,Nt);
    end
    
    % Particle-defined fields:
    if sum(strcmpi('Bx_p',vars.output))
        Bx_p{ss}   = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('dBx_p',vars.output))
        dBx_p{ss}  = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('ddBx_p',vars.output))
       ddBx_p{ss} = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('Ex_p',vars.output))
        Ex_p{ss}   = zeros(N_CP,Nt);
    end
        
end

% Get particle-defined quantities:
% ============================
for rr = 1:ranksParticles
    
    disp(['Particle rank #', num2str(rr),' out of ',num2str(ranksParticles)]);
    
    % File name:
    fileName = ['PARTICLES_FILE_',num2str(rr-1),'.h5'];
    
    % Get data:
    IONS = HDF2Struct_v2(fileName);
    
    % Loop over all ion species:
    for ss = 1:numIonSpecies
        
        % Species field:
        speciesField = ['spp_',num2str(ss)];
              
        % Number of particles per rank:
        N_CP_MPI = double(ionParameters{ss}.NSP_OUT);
                        
        % Data range:
        rng_p = 1 + N_CP_MPI*(rr-1) : 1 : N_CP_MPI*rr;

        % Loop over all time:
        for tt = 1:Nt
            % Obtain time stamp:
            kk = str2double(timeSteps{tt}(2:end)) + 1;
            
            % Time:
            t_p(kk) = IONS.(timeSteps{tt}).time;
            
           % Particle states:
            if sum(strcmpi('x_p',vars.output))
                x_p{ss}(rng_p,kk)    = IONS.(timeSteps{tt}).ions.(speciesField).x_p;
            end
            
            if sum(strcmpi('V_p',vars.output))
                vpar_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).V_p(:,1);
                vper_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).V_p(:,2);   
            end
            
            if sum(strcmpi('a_p',vars.output))
                a_p{ss}(rng_p,kk)    = IONS.(timeSteps{tt}).ions.(speciesField).a_p;
            end
            
            if sum(strcmpi('mu_p',vars.output))
                mu_p{ss}(rng_p,kk)   = IONS.(timeSteps{tt}).ions.(speciesField).mu_p;
            end        
        
            % Particle-defined moments:
            if sum(strcmpi('n_p',vars.output))
                n_p{ss}(rng_p,kk)    = IONS.(timeSteps{tt}).ions.(speciesField).n_p;
            end
            
            if sum(strcmpi('nv_p',vars.output))
                nv_p{ss}(rng_p,kk)   = IONS.(timeSteps{tt}).ions.(speciesField).nv_p;
            end
            
            if sum(strcmpi('Tpar_p',vars.output))
                Tpar_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).Tpar_p;
            end
            
            if sum(strcmpi('Tper_p',vars.output))
                Tper_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).Tper_p;
            end
            
            if sum(strcmpi('Te_p',vars.output))
                Te_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).Te_p;
            end
            
            if sum(strcmpi('u_p',vars.output))
                ux_p{ss}(rng_p,kk)    = IONS.(timeSteps{tt}).ions.(speciesField).u_p.x;               
            end
    
            % Mesh-defined moments:
            if rr == 1
                if sum(strcmpi('n_m',vars.output))
                    n_m{ss}(:,kk)    = IONS.(timeSteps{tt}).ions.spp_1.n_m;
                end
                
                if sum(strcmpi('nv_m',vars.output))
                    nv_m{ss}   = IONS.(timeSteps{tt}).ions.spp_1.nv_m;
                end
                
                if sum(strcmpi('Tpar_m',vars.output))
                    Tpar_m{ss}(:,kk) = IONS.(timeSteps{tt}).ions.spp_1.Tpar_m;
                end
                
                if sum(strcmpi('Tper_m',vars.output))
                    Tper_m{ss}(:,kk) = IONS.(timeSteps{tt}).ions.spp_1.Tper_m;
                end
                
                if sum(strcmpi('Te_m',vars.output))
                    Te_m{ss}(:,kk) = IONS.(timeSteps{tt}).ions.spp_1.Te_m;
                end
                
                if sum(strcmpi('u_m',vars.output))
                    ux_m{ss}(:,kk)   = IONS.(timeSteps{tt}).ions.spp_1.u_m.x;                
                end
            end
                
            % Particle-defined fields:
            if sum(strcmpi('Bx_p',vars.output))
                Bx_p{ss}(rng_p,kk)   = IONS.(timeSteps{tt}).ions.(speciesField).Bx_p;
            end
            
            if sum(strcmpi('dBx_p',vars.output))
                dBx_p{ss}(rng_p,kk)  = IONS.(timeSteps{tt}).ions.(speciesField).dBx_p;
            end
            
            if sum(strcmpi('ddBx_p',vars.output))
                ddBx_p{ss}(rng_p,kk) = IONS.(timeSteps{tt}).ions.(speciesField).ddBx_p;
            end
            
            if sum(strcmpi('Ex_p',vars.output))
                Ex_p{ss}(rng_p,kk)   = IONS.(timeSteps{tt}).ions.(speciesField).Ex_p;
            end         
                                    
        end
        
    end
end
cd ..

disp('Particle data extraction completed!')
disp('********************************************************************')
disp('')

%% Testing the RK4 integrator results:
if 0
    close all
    rng = find(x_p{1}(:,1) >1.50 & x_p{1}(:,1)< 1.55);
%     rng = find(x_p{1}(:,1) >0.0 & x_p{1}(:,1)< 3.55);

    % For symmetric domains:
     rng = find(x_p{1}(:,1) >-0.1 & x_p{1}(:,1)< 0.1);


    % Kinetic energy:
    Ma   = m.ions.spp_1.M;
    KE_p = 0.5*Ma*(vpar_p{1}.^2 + vper_p{1}.^2)/e_c;
    
    % Magnetic moment:
    mu_pb = (Ma*vper_p{1}.^2)./(2*Bx_p{1});
    
    % Time:
    tMax = max(t_p);
    
    % Magnetic field:
    bMax = max(Bx_m(:,1));
    
    figure;
    hold on
    plot(x_m,Bx_m(:,1)*tMax/bMax)
    for ii = 1:numel(rng)
        plot(x_p{1}(rng(ii),:),t_p,'k-'); 
    end
    xlim([0,m.geometry.LX]);
    title('Particle position vs time')
    
    figure;
    hold on
    for ii = 1:numel(rng)
        plot(t_p,KE_p(rng(ii),:)); 
    end
    xlim([0,tMax]);
    title('Particle KE vs time')
    
    try
        figure;
        hold on
        for ii = 1:numel(rng)
            plot(t_p,mu_p{1}(rng(ii),:)); 
        end
        xlim([0,tMax]);
        title('mu vs time')    
    end
    
    figure; 
    hold on
    for ii = 1:numel(rng)
        plot3(x_p{1}(rng(ii),:),t_p,KE_p(rng(ii),:)); 
    end
    xlim([0,m.geometry.LX]);
    zlim([1,1e2])
    title('Kinetic energy vs position')

    figure; 
    hold on
    for ii = 1:numel(rng)
        plot3(x_p{1}(rng(ii),:),t_p,mu_pb(rng(ii),:)); 
    end
    xlim([0,m.geometry.LX]);
    for ii = 1:numel(rng)
        plot3(x_p{1}(rng(ii),:),t_p,mu_p{1}(rng(ii),:),'.'); 
    end
    xlim([0,m.geometry.LX]);
    title('Magnetic moment vs position')  
    
    % Test individual particles:
    figure; 
    ii = 10; 
    subplot(1,3,1); 
    plot(t_p,x_p{1}(ii,:),'r.-'); 
    ylim([0,3])
    
    subplot(1,3,2); 
    plot(t_p,KE_p(ii,:)); 
    ylim([0,1e3]); 
    
    subplot(1,3,3); 
    plot(t_p,mu_p{1}(ii,:)); 
    ylim([0,1e-16])
    xlim([0,t_p(end)])
end

%% Plot ion moments:
close all

if find(strcmpi('Bx_p',vars.output) == 1) && find(strcmpi('x_p',vars.output) == 1)
    figure
    hold on
    plot(x_p{1},Bx_p{1},'k.')
    plot(x_m   ,Bx_m,'r.-')
end

if find(strcmpi('n_m',vars.output) == 1)
    figure
    mesh(t_p,x_m,movmean(n_m{1},10,1));
    zlim([0,2e20])
    caxis([0,2e20])
    title('n')
end

if find(strcmpi('Ex_m',vars.output) == 1)
    figure
    Te = double(m.Te);
    mesh(t_p,x_m,movmean(movmean(Ex_m/Te,10,1),4,2));
    zlim([-1,1]*8)
    caxis([-1,1]*8)
    title('Ex')
end

if find(strcmpi('u_m',vars.output) == 1)
    Cs = sqrt( e_c*(Te_m{1} + 3*Tpar_m{1})/m.ions.spp_1.M);
    figure
    mesh(t_p,x_m,movmean(ux_m{1}./Cs,10,1));
    title('u_x/{v_T}')
end

if find(strcmpi('Tpar_m',vars.output) == 1)
    figure
    mesh(t_p,x_m,movmean(Tpar_m{1},10,1));
    Tpar_max = max(max(Tpar_m{1}));
    zlim([0,2*Tpar_max])
    title('T_par')
end

if find(strcmpi('Tper_m',vars.output) == 1)
    figure
    mesh(t_p,x_m,movmean(Tper_m{1},10,1));
    zlim([0,1000])
    title('T_per')
end

% Cross sectional area:
A0 = pi*(5/100).^2;
B0 = 0.2;
A  = A0*B0./Bx_m;

if find(strcmpi('u_m',vars.output) == 1)
    % Plasma integrated flux:
    F = n_m{1}.*ux_m{1}.*A;

    figure
    mesh(t_p,x_m,movmean(F,10,1));
    title('particle flux')
end

% Temperature during and before RF
figure('color','w')
hold on
rng = 18:20; 
hT(1) = plot(x_m,movmean(mean(Tpar_m{1}(:,rng),2),9,1)); 
hT(2) = plot(x_m,movmean(mean(Tper_m{1}(:,rng),2),9,1));
rng = 30:40; 
hT(3) = plot(x_m,movmean(mean(Tpar_m{1}(:,rng),2),9,1)); 
hT(4) = plot(x_m,movmean(mean(Tper_m{1}(:,rng),2),9,1));
plot(x_m,Bx_m(:,1)*200)


%% Force balance:
% close all

% x range:
rng_x = find(x_m >= -2 & x_m <= 10);

% Plasma quantities:
M      = m.ions.spp_1.M;
A0     = pi*(5/100)^2;
B0     = 0.2;
Te     = Te_m{1}(rng_x,:);
Ti_par = Tpar_m{1}(rng_x,:);
Ti_per = Tper_m{1}(rng_x,:);
ne     = n_m{1}(rng_x,:);
U      = ux_m{1}(rng_x,:);
B      = Bx_m(rng_x,:);
A      = A0*B0./B;
x      = x_m(rng_x,:);
t      = t_p;
dx     = diff(x(1:2));
dt     = diff(t(1:2));

% Derived:
flux_density  = ne.*U;
flux          = flux_density.*A;
E_mean        = e_c*(Te + Ti_par + Ti_per);
power_density = flux_density.*E_mean;
power         = power_density.*A;
powerFlux_par = e_c*(Te + Ti_par + Ti_per + 13.6).*flux_density;

% Smooth in space:
fr = 3;
ne     = movmean(ne    ,fr,1);
Ti_par = movmean(Ti_par,fr,1);
Ti_per = movmean(Ti_per,fr,1);
U      = movmean(U     ,fr,1);

% Array size:
Nx = numel(x);
Nt = numel(t);

% Plasma flux density:
nU = ne.*U;

% Pressure:
P_par = e_c*(Te + Ti_par).*ne;
P_per = e_c*(Te + Ti_per).*ne;
P_KE  = M*ne.*U.^2; 

% Smooth over space:
fr = 15;
P_par = movmean(P_par,fr,1);
P_per = movmean(P_per,fr,1);
P_KE  = movmean(P_KE ,fr,1);
nU    = movmean(nU'  ,fr,1)';

% Gradients (central difference):
dP_par_dx = central_diff(P_par)/dx;
dB_dx     = central_diff(B)/dx;
dnu_dt    = (central_diff(nU')/dt)';

% Derived:
dP = P_per - P_par;

% Forces:
F_KE_x = B.*central_diff(P_KE./B)/dx;
F_KE_t = M*dnu_dt;
F_par  = -dP_par_dx;
F_mag  = -(dP./B).*dB_dx;

% Smooth over time:
fr = 5;
F_KE_x = movmean(F_KE_x,fr,2);
F_KE_t = movmean(F_KE_t,fr,2);
F_par  = movmean(F_par ,fr,2);
F_mag  = movmean(F_mag ,fr,2);

figure; 
rng = 26:41;
% rng = 50:60;
plot(x,mean(P_KE(:,rng),2))

figure; 
plot(x,mean(dP_par_dx(:,rng),2))

figure('color','w');
hold on
% rng = (Nt-5): Nt;
hF(1) = plot(x,mean(F_KE_x(:,rng),2),'k','LineWidth',2);
hF(2) = plot(x,mean(F_par(:,rng) ,2),'bl','LineWidth',2);
hF(3) = plot(x,mean(F_mag(:,rng) ,2),'r','LineWidth',2);
box on
hL = legend(hF,'$f_{K}$','$f_{\parallel}$','$f_{B}$');
set(hL,'interpreter','latex','fontSize',14)
ylim(3*[-1,1]*max(mean(F_KE_x(:,rng),2)))
xlabel('x [m]','interpreter','latex','fontSize',14)
ylabel('[Nm$^{-3}$]','interpreter','latex','fontSize',14)
grid on
title(['ICRF ',myDir(6:(end-7)),' [kW]'],'interpreter','latex','fontSize',14)

figure('color','w'); 
hold on
box on
f = F_KE_t(:,rng) + F_KE_x(:,rng);
g = F_par(:,rng) + F_mag(:,rng);
hL(1) = plot(x,+mean(f,2),'k','LineWidth',2)
hL(2) = plot(x,-mean(g,2),'r','LineWidth',2)
ylim(3*[-1,1]*max(mean(F_KE_x(:,rng),2)))
title(['Paralle force balance'],'interpreter','latex','fontSize',14)
hLeg = legend(hL,'$f_{K}$','$f_{\parallel} + f_{B}$');
set(hLeg,'interpreter','latex','fontSize',14)

%% Plots to compare with SOLPS

% close all

% Plasma density:
figure('color','w')
hold on
box on
g = n_m{1}(:,end-30:end);
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,1e20])
title('Plasma density [m$^{-3}]$','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)


% Electron tempeture:
figure('color','w')
hold on
box on
g = Te_m{1}(:,end-10:end);
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,6])
title('Electron temperature [eV]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Plasma flow:
figure('color','w')
hold on
box on
Cs = sqrt( e_c*(Te_m{1} + 3*Tpar_m{1})/m.ions.spp_1.M);
f = ux_m{1}./Cs;
Mach = mean(f(:,end-30:end),2);
fmax = max(Mach);
Bmax = max(Bx_m(:,1));
plot(x_m,Mach,'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([-2.5,2])
title('Plasma parallel flow $U_\parallel/v_T$','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Electron pressure:
figure('color','w')
hold on
box on
g = n_m{1}(:,end-30:end).*Te_m{1}(:,1)*e_c;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,1]*50)
title('Electron pressure [Pa]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Parallel energy flux:
figure('color','w')
hold on
box on
g = powerFlux_par*1e-6;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([-3,3])
title('Parallel Energy flux [MWm$^{-2}$]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Total power flow:
figure('color','w')
hold on
box on
g = powerFlux_par.*A;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([-1,1]*6e3)
title('Parallel power flow [W]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Convective heat flux:
figure('color','w')
hold on
box on
g = flux_density.*(Ti_par + Ti_per)*e_c*1e-6;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([-1,1]*0.8)
title('Parallel convective heat flux [MWm$^{-2}$]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

return

%% Calculate velocity space PDF:
% Select range:

rangeType = 1;
switch rangeType
    case 1
        % mirror region:
        % % ===============
        x_center = 1.6;
        x_delta  = 0.8;
        z1 = x_center - x_delta;
        z2 = x_center + x_delta;
    case 2
        % exhaust region:
        % ===============
        z1 = 1.6;
        z2 = 2.5;
    case 3
        % MPEX RF region:
        % ===============
        z1 = 5.1;
        z2 = 5.6;   
    case 4
        % MPEX Target region:
        % ===============
        z1 = 7.0;
        z2 = 7.9;   
    case 5
        % MPEX source region:
        % ===============
        z1 = -0.25 + 1;
        z2 = 0.25 + 1;
end

% Ion parameters:
beam.E = 1.5e3;
Ma = m.ions.spp_1.M;

% Ion derived quantities:
vT = double(sqrt(m.Te*e_c/Ma));

% Create mesh with ghost cells"
Nx = double(m.geometry.NX_IN_SIM);
dx = double(m.geometry.DX);
x_m_g = ((1:(Nx+4))-1)*dx + 0.5*dx - 2*dx;

% Select range of velocity grid:
a = 12;
vxMin = -a*vT;
vxMax = +a*vT;
vyMin = -a*vT;
vyMax = +a*vT;

NV = 200;
dvx = (vxMax-vxMin)/NV;
dvy = (vyMax-vyMin)/NV;
vxGrid = ((1:NV)-1)*dvx + dvx/2 + vxMin;
vyGrid = ((1:NV)-1)*dvy + dvy/2 + vyMin;

% Number of time steps:
NS = size(a_p{1},2);

% Allocate memory:
fv = zeros(NV+4,NV+4,NS);

% Create magnetic field vector on mesh with ghost cells:
% Need to have valid magnetic field data in the ghost cells
% B = zeros(numel(x_m_g), 1);
% B = interp1(x_m,Bx_m(:,1),x_m_g,'spline','extrap')';
B0 = Bx_m(round(numel(Bx_m(:,1))/2),1);
A0 = pi*(0.1^2);
% A = B0*A0./B;

Phi0 = B0*A0;
for jj = 1:NS
    % Select the spatial range:
    rng = find(x_p{1}(:,jj) > z1 & x_p{1}(:,jj) < z2);
%     rng = rng(1:150E3);
    Nrng = numel(rng);
    Rm = rand(Nrng,1);
        
    vvy = vper_p{1}(rng,jj).*cos(2*pi*Rm);
    vvx = vpar_p{1}(rng,jj);   
    
    [WXL,WXC,WXR,MXI,WYL,WYC,WYR,MYI] = AssignCell_2D(vvx/vT,vvy/vT,vxGrid/vT,vyGrid/vT);
    disp(['jj: ',num2str(jj)])
    
    for ii = 1:Nrng
        if ~isempty(MXI(ii))||~isempty(MYI(ii))
            ix = MXI(ii) + 2;
            iy = MYI(ii) + 2;

            fv(ix - 1,iy + 0,jj) = fv(ix - 1,iy + 0,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXL(ii)*WYC(ii)*a_p{1}(rng(ii),jj)/Nrng;
            fv(ix + 0,iy + 0,jj) = fv(ix + 0,iy + 0,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXC(ii)*WYC(ii)*a_p{1}(rng(ii),jj)/Nrng;    
            fv(ix + 1,iy + 0,jj) = fv(ix + 1,iy + 0,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXR(ii)*WYC(ii)*a_p{1}(rng(ii),jj)/Nrng;

            fv(ix + 0,iy - 1,jj) = fv(ix + 0,iy - 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXC(ii)*WYL(ii)*a_p{1}(rng(ii),jj)/Nrng;
            fv(ix + 0,iy + 1,jj) = fv(ix + 0,iy + 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXC(ii)*WYR(ii)*a_p{1}(rng(ii),jj)/Nrng;    

            fv(ix + 1,iy - 1,jj) = fv(ix + 1,iy - 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXR(ii)*WYL(ii)*a_p{1}(rng(ii),jj)/Nrng;
            fv(ix + 1,iy + 1,jj) = fv(ix + 1,iy + 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXR(ii)*WYR(ii)*a_p{1}(rng(ii),jj)/Nrng;

            fv(ix - 1,iy - 1,jj) = fv(ix - 1,iy - 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXL(ii)*WYL(ii)*a_p{1}(rng(ii),jj)/Nrng;
            fv(ix - 1,iy + 1,jj) = fv(ix - 1,iy + 1,jj) + (Bx_p{1}(rng(ii),jj)/Phi0)*WXL(ii)*WYR(ii)*a_p{1}(rng(ii),jj)/Nrng;
        end
    end
end


%% Plot PDF:
% Mirror ratio
R_m =max(Bx_m(:,1))/B0;

% Loss cone angle:
theta_m = asin(1/sqrt(R_m));
theta_NBI = 45*pi/180;
vNBI = sqrt(e_c*5000/Ma);

% % Plotting the integral of the electric field:
figure; 
endTime = size(Tpar_m{1},2);
rngLC = endTime-10:endTime-1;
frame=1;
ff = mean(movmean(Ex_m(:,rngLC),frame,1),2);
try
    rng = find(x_m>x_center & x_m<z2);
catch
        rng = find(x_m>z1 & x_m<z2);
end

delta_V = trapz(x_m(rng),ff(rng));
gg = Bx_m(:,jj);
plot(x_m(:),ff,'k');hold on;
plot(x_m(rng),ff(rng),'g');
plot(x_m(:),gg*max(ff)/max(gg),'r');
ylim([-20 20]*m.Te);
hold off;

% Calculating the hyperbola that defines the trapped region:
vvpar = vxGrid/vT;
vvper = sqrt((vvpar.^2 + 1*delta_V/(m.Te))/(R_m - 1));

figure('color','w'); 
rngx = 3:(NV+4-2);
rngy = 3:(NV+4-2);
Limit=10;
fmax = max(max(fv(:,:,end)));
for jj = 1:1:size(fv,3)-5
    contourf(vxGrid/vT,vyGrid/vT,mean(fv(rngx,rngy,jj:jj+5),3)',10,'LineStyle','none');

    hold on
    lossCone(1) = line([-1,+1]*Limit,[-1,+1]*Limit*tan(theta_m));
    lossCone(2) = line([-1,+1]*Limit,[+1,-1]*Limit*tan(theta_m));
    set(lossCone,'LineWidth',2,'color','k','lineStyle','--','LineWidth',2)
        
    if 0
        NBICone(1) = line([-1,+1]*Limit,[-1,+1]*Limit*tan(theta_NBI));
        NBICone(2) = line([-1,+1]*Limit,[+1,-1]*Limit*tan(theta_NBI));
        set(NBICone,'LineWidth',2,'color','g','lineStyle','--','LineWidth',2)
        plot(vNBI*cos(theta_NBI)/vT,vNBI*sin(theta_NBI)/vT,'go')
    end
    
    
    confinementRegion(1) = plot(vvpar,+vvper,'k','LineWidth',1);
    confinementRegion(2) = plot(vvpar,-vvper,'k','LineWidth',1);
    
    title(['$f(v_\parallel,v_\perp)$ '],'interpreter','latex','FontSize',20)
    hold off
    axis square
    xlim([-1,1]*Limit)
    ylim([-1,1]*Limit)
    xlabel('$v_\parallel/v_{Ti}$','interpreter','latex','FontSize',20)
    ylabel('$v_\perp/v_{Te}$','interpreter','latex','FontSize',20)
    caxis([0,0.9]*fmax)
    colorbar
    colormap(flipud(hot))
    view([0,90])
    drawnow
    pause(0.01)
end


return

%%
close all
figure
hold on
rng = 136;
phi = -cumtrapz(x_m,ff);
phi = phi - max(phi);
for ii = rng
    plot(t_p,x_p{1}(ii,:))
    KE{ii} = 0.5*m_p*( vper_p{1}(ii,:).^2 + vpar_p{1}(ii,:).^2 )/e_c;
    KE_per{ii} = 0.5*m_p*( vper_p{1}(ii,:).^2 )/e_c;
    KE_par{ii} = 0.5*m_p*( vpar_p{1}(ii,:).^2 )/e_c;
    phi_p{ii}  = double(interp1(x_m,phi,x_p{1}(ii,:)));
    mm_p{ii} = 0.5*m_p*(vper_p{1}(ii,:).^2)./Bx_p{1}(ii,:);
end

figure
hold on
for ii = rng
plot3(t_p,x_p{1}(ii,:),KE{ii},'g.-')
plot3(t_p,x_p{1}(ii,:),-phi_p{ii},'m.-')

plot3(t_p,x_p{1}(ii,:),KE{ii} + phi_p{ii},'k.-')
end
zlim([0,30])
ylim([-2,8])

figure
hold on
for ii = rng
    plot3(t_p,x_p{1}(ii,:),KE_par{ii} ,'k.-')
    plot3(t_p,x_p{1}(ii,:),KE_per{ii} ,'r.-')
end
zlim([0,30])
ylim([-2,8])

figure
hold on
for ii = rng
    plot3(t_p,x_p{1}(ii,:),mm_p{ii},'k.-')
end
ylim([-2,8])
zlim([0,1.3*max(mm_p{ii})])
%%
% Animation of the density change:
figure('color','w')
for jj = 1:numel(t_p)
    plot(x_m,movmean(n_m{1}(:,jj),10,1))
    hold on
    plot(x_m,Bx_m(:,1)*(1e20)/3)
    hold off
    ylim([0,1.2e20]);
    title(['frame = ',num2str(jj)])
    drawnow
    pause(0.1)
end
hold on
plot(x_m,movmean(n_m{1}(:,20),10,1))
plot(x_m,movmean(n_m{1}(:,26),10,1))
plot(x_m,movmean(n_m{1}(:,30),10,1))
% plot(x_m,movmean(n_m{1}(:,40),10,1))

% Animation of the flow change:
figure('color','w')
for jj = 1:numel(t_p)
    plot(x_m,movmean(ux_m{1}(:,jj)./Cs(:,jj),10,1))
    hold on
    plot(x_m,Bx_m(:,1))
    hold off
    ylim([-2.5,2.5]);
    title(['frame = ',num2str(jj)])
    grid on
    drawnow
    pause(0.1)
end
hold on
plot(x_m,movmean(ux_m{1}(:,20)./Cs(:,20),10,1))
plot(x_m,movmean(ux_m{1}(:,26)./Cs(:,20),10,1))
plot(x_m,movmean(ux_m{1}(:,30)./Cs(:,20),10,1))

% Animation of the temperature change:
figure('color','w')
for jj = 1:numel(t_p)
    hT(1) = plot(x_m,movmean(Tpar_m{1}(:,jj),10,1));
    hold on
    hT(2) = plot(x_m,movmean(Tper_m{1}(:,jj),10,1));
    plot(x_m,Bx_m(:,1)*10);
    hold off
    ylim([0,20]);
    title(['frame = ',num2str(jj)])
    grid on
    drawnow
    pause(0.1)
end
legend(hT,'T_{\perp}','T_{\perp}')


%% Calculate moments:

% Calculate mesh-defined ion moments:
% =========================================================================
Nx_pdf = 200;

% Loop over all species:
for ss = 1:numIonSpecies    
    
    % Ion parameters:
    M = ionParameters{ss}.M;
    
    % Real particles for each computational particle:
    K = ionParameters{ss}.K;

    % Total number of computational particles:
    N_CP = ionParameters{ss}.N_CP*ranksParticles;

    % Total number of real particles represented:
    NR = N_CP*K;
        
    % Characteristic temperatures:
    Tpar = ionParameters{ss}.Tpar;
    Tper = ionParameters{ss}.Tper;

    % Loop over all time:
    for tt = 1:Nt
        
        % "x" grid:
        x_max = double(m.geometry.LX);
        xbin = linspace(0,x_max,Nx_pdf);

        % "vx" grid:
        vx_max = sqrt(1*e_c*Tpar/M);
        vxbin = linspace(-12,+12,Nx_pdf)*vx_max;

        % "vy" grid:
        vy_max = sqrt(1*e_c*Tper/M);
        vybin = linspace(-12,+12,Nx_pdf)*vy_max;

        % "vz" grid:
        vz_max = vy_max;
        vzbin = vybin;    

        % Inputs for PDF calculation:
        x  = xbin(1:end-1);
        vx = vxbin(1:end-1);
        vy = vybin(1:end-1);
        vz = vybin(1:end-1);

        % Calculate compression factor:
        Bx = interp1(x_m,Bx_m.x(:,tt),x)'; 
        % Select model:
        modelType = 2;
        switch modelType
            case 1
                compFactor = ones(size(Bx))';
            case 2
                B0 = Bx(round(length(x)/2));
                compFactor = double(Bx/B0);
        end
        [~,CP] = meshgrid(compFactor(:,1));

        % Calculate PDFs:
        fxvx{ss}{tt} = CP.*histcounts2(x_p{ss}(:,tt),vx_p{ss}(:,tt),xbin,vxbin,'Normalization','pdf');
        fxvy{ss}{tt} = CP.*histcounts2(x_p{ss}(:,tt),vy_p{ss}(:,tt),xbin,vybin,'Normalization','pdf');
        fxvz{ss}{tt} = CP.*histcounts2(x_p{ss}(:,tt),vz_p{ss}(:,tt),xbin,vzbin,'Normalization','pdf');

        rng=find(x_p{ss}(:,tt)>1 & x_p{ss}(:,tt)<2.2);
%         vper = sqrt(vy_p{ss}(rng,tt).^2 + vz_p{ss}(rng,tt).^2);
        vper = vz_p{ss}(rng,tt);
        fvxvper{ss}{tt} = histcounts2(vx_p{ss}(rng,tt),vper,vxbin,vzbin,'Normalization','pdf');

        % Zeroth moment:
        n{ss}(:,tt)  = NR*trapz(vx,fxvx{ss}{tt}.*vx.^0,2);
        
        % First moment:
        ux{ss}(:,tt) = NR*trapz(vx,fxvx{ss}{tt}.*vx.^1,2)./n{ss}(:,tt);
        uy{ss}(:,tt) = NR*trapz(vy,fxvy{ss}{tt}.*vy.^1,2)./n{ss}(:,tt);
        uz{ss}(:,tt) = NR*trapz(vz,fxvz{ss}{tt}.*vz.^1,2)./n{ss}(:,tt);
        
        % Second moment:
        Px{ss}(:,tt) = M*NR*trapz(vx,fxvx{ss}{tt}.*(vx - ux{ss}(:,tt)).^2,2);
        Py{ss}(:,tt) = M*NR*trapz(vy,fxvy{ss}{tt}.*(vy - uy{ss}(:,tt)).^2,2);
        Pz{ss}(:,tt) = M*NR*trapz(vz,fxvz{ss}{tt}.*(vz - uz{ss}(:,tt)).^2,2);
        
        % Temperature:
        Tx{ss}(:,tt) = Px{ss}(:,tt)./(n{ss}(:,tt)*e_c);
        Ty{ss}(:,tt) = Py{ss}(:,tt)./(n{ss}(:,tt)*e_c);
        Tz{ss}(:,tt) = Pz{ss}(:,tt)./(n{ss}(:,tt)*e_c);
        
    end
end

% Calculate mesh-defined electron moments:
% =========================================================================
% Initialize variables:
ne_m  = zeros(Nx,1);
nue_m = zeros(Nx,1);

% Calculate electron density and flux density:
for ss = 1:numIonSpecies    
    Z = ionParameters{ss}.Z;
    ne_m = ne_m + n_m{ss}*Z;
    nue_m = nue_m + n_m{ss}.*ux_m{ss}*Z;
end

% Calculte electron drift velocity and pressure:
Ux_m = nue_m./ne_m;
Pe_m = ne_m*m.Te*k_B/e_c;

%% Plot data:
% =========================================================================
close all
 
% Select ion species:
ss = 1;

% Select time:
% jj =249;
% jj = 5;
% jj = 10;
% jj = 50;
jj = 11;
% jj = 60;
% jj = 10;

% Plot distribution function:
% =========================================================================
figure('color','w');
subplot(2,2,1)
surf(xbin(1:end-1),vxbin(1:end-1)./vx_max,fxvx{ss}{jj}','LineStyle','none');
title('x-v_x');
xlabel('x')
ylabel('v_x')
view([0 90])
xlim([0 3.2]);
ylim([-6 6]);

subplot(2,2,2)
surf(xbin(1:end-1),vybin(1:end-1)./vy_max,fxvy{ss}{jj}','LineStyle','none');
title('x-v_y')
xlabel('x')
ylabel('v_y')
view([0 90])
xlim([0 3.2]);
ylim([-6 6]);

subplot(2,2,3)
surf(xbin(1:end-1),vzbin(1:end-1)./vy_max,fxvz{ss}{jj}','LineStyle','none');
title('x-v_z');
xlabel('x')
ylabel('v_z')
view([0 90])
xlim([0 3.2]);
ylim([-6 6]);

subplot(2,2,4)
surf(vxbin(1:end-1)./vx_max,vzbin(1:end-1)./vz_max,fvxvper{ss}{jj}','LineStyle','none');
title('v_x-v_{per}');
xlabel('v_x')
ylabel('v_{per}')
view([0 90])
xlim([-4 4]*3);
ylim([-4 4]*3);

set(gcf,'position',[295   597   560   420]);

% Ex, Ux, U parallel and magnetic field:
% =========================================================================
% Select time:
% rngTime = 15:20;
% rngTime = 10;
rngTime = 9:11;

% Smooth data:
Ex0     = mean(E_m.x(:,rngTime),2);
Ex_mean = movmean(Ex0,40);

% Calculate error bars:
for ii = 1:numel(rngTime)
    dEx(:,ii)  = movmean(E_m.x(:,rngTime(ii)),2) - Ex_mean;
end
Ex_std = 0.5*movmean(std(dEx,1,2),40);

% Ex:
% -------------------------------------------------------------------------
figure('color','w')
subplot(2,1,1)
hold on
Te = m.Te*k_B/e_c;
plot(x_m,Ex_mean/Te,'k','LineWidth',2)
plot(x_m,(Ex_mean + Ex_std)/Te,'k:','LineWidth',1.5)
plot(x_m,(Ex_mean - Ex_std)/Te,'k:','LineWidth',1.5)
set(gca,'FontName','times','fontSize',18)
ylabel('$E_\parallel / T_e$','interpreter','latex','fontSize',22)
xlabel('$x$ [m]','interpreter','latex','fontSize',22)
box on
ylim([-1,1]*20)
xlim([0,max(x_m)])

% Ux:
% -------------------------------------------------------------------------

% Select time:
% rngTime = 20;
rngTime = 10;
% rngTime = 13;

subplot(2,1,2)
hold on
yyaxis left
ax(1) = gca;
hb(2) = plot(x_m,Ux_m(:,rngTime)/vx_max,'k','LineWidth',2);
set(gca,'YColor','k','FontName','times','fontSize',18)
ylabel(ax(1),'$u_\parallel / u_{T_\parallel}$','interpreter','latex','fontSize',22)
xlabel('$x$ [m]','interpreter','latex','fontSize',22)

yyaxis right
ax(2) = gca;
hb(1) = plot(x_m,Bx_m.x(:,rngTime),'r','LineWidth',2');
ylabel(ax(2),'$B_0$','interpreter','latex','fontSize',22)
set(ax(2),'YColor','k')
box on
% Formatting:
ylim(ax(1),[-10,10])
ylim(ax(2),[-1,1]*20)
xlim(ax,[0,max(x_m)])
hL = legend(hb,'$B_0$','$u_\parallel / u_{T_\parallel}$');
set(hL,'interpreter','latex','location','southeast','FontSize',18)

set(gcf,'position',[934   597   560   420]);

% Ion Temperatures:
% =========================================================================
figure;
subplot(3,1,1)
hold on
frame=30;
hT(1) = plot(xbin(1:end-1),movmean(Tx{ss}(:,jj),frame),'k');
hT(2) = plot(xbin(1:end-1),movmean(Ty{ss}(:,jj),frame),'r');
hT(3) = plot(xbin(1:end-1),movmean(Tz{ss}(:,jj),frame),'g');
legend(hT,'Tx','Ty','Tz')
title('Temperature, t= $3000 \Omega_c^{-1}$','interpreter','latex','fontSize',18)
ylim([0,2500])
box on;

ylabel('$T_i$ [eV]','interpreter','Latex','fontSize',14)
xlim([0 3.2]);
grid on

% Ion and electron density:
% =========================================================================
subplot(3,1,2)
hold on
plot(x_m(1:length(ne_m(:,jj))),ne_m(:,jj),'k');
plot(x_m(1:length(ne_m(:,jj))),n_m{ss}(:,jj),'r');
title('density, t= $3000 \Omega_c^{-1}$','interpreter','latex','fontSize',18)
box on;

ylabel('$n_i$ [m$^{-3}$]','interpreter','Latex','fontSize',14)
xlim([0 3.2]);
ylim([0 6E19])
grid on

subplot(3,1,3)
plot(x_m(1:length(Bx_m.x(:,jj))),Bx_m.x(:,jj),'b');hold on;
plot(x_m(1:length(Bx_m.y(:,jj))),Bx_m.y(:,jj),'r')
title('Magnetic field, t= $3000 \Omega_c^{-1}$','interpreter','latex','fontSize',18)
xlabel('$x$ [$m$]','interpreter','Latex','fontSize',14)
ylabel('$B$ [T]','interpreter','Latex','fontSize',14)
hL = legend('Axial B-field','Radial B-field');
set(hL,'interpreter','latex','location','southeast','FontSize',14)
grid on
xlim([0 3.2]);
set(gcf,'position',[ 406    44   560   420]);

%% Functions:

function f = central_diff(y)
    Nx = size(y,1);
    f = zeros(size(y));
    f(2:(Nx-1),:) = 0.5*( y(3:Nx,:) - y(1:(Nx-2),:) );
    f(1,:)  = f(2,:);
    f(Nx,:) = f(Nx-1,:);
end

function [WXL,WXC,WXR,MXI,WYL,WYC,WYR,MYI] = AssignCell_2D(Xi,Yi,Xgrid,Ygrid)
% Xgrid and Ygrid are uniform grids

% Total number of computational particles:
NC = numel(Xi);

% Size of grids:
Nx = numel(Xgrid);
NY = numel(Ygrid);

% Grid spacing:
dX = Xgrid(2)-Xgrid(1);
dY = Ygrid(2)-Ygrid(1);

% Grid offset:
Xoffset = min(Xgrid) - 0.5*dX;
Yoffset = min(Ygrid) - 0.5*dY;

% Allocate memory:
MXI = zeros(size(Xi));
WXL = zeros(size(Xi));
WXC = zeros(size(Xi));
WXR = zeros(size(Xi));

MYI = zeros(size(Xi));
WYL = zeros(size(Xi));
WYC = zeros(size(Xi));
WYR = zeros(size(Xi));


for ii = 1:NC
    % Find the index "MXI" of nearest grid point to Xi:
    xii = Xi(ii) - Xoffset;
    nnx = xii/dX;
    MXI(ii) = round(nnx + 1/2);
    if MXI(ii) <= Nx && MXI(ii) > 1
        X = Xgrid(MXI(ii)) - Xi(ii);
    else 
        MXI(ii) = [];
    end
   
    try
    if ~isempty(MXI(ii))
        WXC(ii) = 0.75 - (X/dX)^2;
        WXL(ii) = 0.5*(1.5 + ((X - dX)/dX) )^2;
        WXR(ii) = 0.5*(1.5 - ((X + dX)/dX) )^2;     
    end
    catch
        disp('hello X')
    end
    
    % Find the index "MYI" of nearest grid point to Yi:
    yii = Yi(ii) - Yoffset;
    nny = yii/dY;
    MYI(ii) = round(nny + 1/2);
    if MYI(ii) <= NY && MYI(ii) > 1
        Y = Ygrid(MYI(ii)) - Yi(ii);
    else 
        MYI(ii) = [];
    end
try    
    if ~isempty(MYI(ii))
        WYC(ii) = 0.75 - (Y/dY)^2;
        WYL(ii) = 0.5*(1.5 + ((Y - dY)/dY) )^2;
        WYR(ii) = 0.5*(1.5 - ((Y + dY)/dY) )^2;        
    end
catch
    disp('hello')
end
end
end

function data=HDF2Struct_v2(f,verbose)
%HDF2STRUCT - Reads HDF5 files into structure
%
% Syntax:  data = HDF2Struct(file_name, verbose)
%
% Inputs:
%    file_name - String. Path to the hdf5 file
%    verbose   - Boolean. Whether or not to print warnings when renaming
%    variables with invalid matlab names
%
% Outputs:
%    data - Matlab structure containing the read data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Luca Amerio
% email: lucamerio89@gmail.com
% March 2019; Last revision: 22-March-2019
%------------- BEGIN CODE --------------
if nargin<2
    verbose=false;
end
data=struct;
loadcontent('/');
    function loadcontent(pathStr)
        % Gets info of current group (or root)
        info=h5info(f,pathStr);
        
        % Loading variables
        for vari=1:length(info.Datasets)
            var=info.Datasets(vari).Name;
            fields  = strsplit(pathStr,'/');   % Ignore initial blank field later
            fields(cellfun(@isempty,fields))=[];
            
            % Validate variable name:
            varName=validateFieldName(var);
            
            % Validate field name:
            fieldsName=validateFieldName(fields);
            
            % Assign data to "data" structure:
            data = setfield(data,fieldsName{:},varName{:},h5read(f,[pathStr,'/',var]));
        end
        
        % Loading attributes
        for atti=1:length(info.Attributes)
            att=info.Attributes(atti).Name;
            fields  = strsplit(pathStr,'/');   % Ignore initial blank field later
            fields(cellfun(@isempty,fields))=[];
            attName=validateFieldName(att);
            fieldsName=validateFieldName(fields);
            data = setfield(data,fieldsName{:},attName{:},h5readatt(f,pathStr,att));
        end
        
        %Loading groups (calls loadcontent recursively for the selected
        %group)
        for gi=1:length(info.Groups)
            loadcontent(info.Groups(gi).Name);
        end
        
        % HDF naming convention allows names unsupported by matlab. This 
        % funtion tryies to clean them when pox_mible.
        function name=validateFieldName(name)
            if ischar(name)
                name={name};
            elseif iscellstr(name)
            else
                error('Input must be either a string or a cell array of strings')
            end
            
            check=~cellfun(@isvarname,name);
            
            if any(check)
                if verbose
                    warning('"%s" is not a valid field name\n',name{check})
                end
                for i=find(check)
                    if any(name{i}==' ')
                        name_new=strrep(name{i},' ','');
                        if verbose
                            warning('"%s" is not a valid field name\nchanging "%s" to "%s"\n',name{i},name{i},name_new)
                        end
                        name{i}=name_new;
                    elseif isnumeric(str2num(name{i}))
                        name{i} = ['t',name{i}];
                    else
                        error('"%s" is not a valid field name\n',name{i})
                    end
                end
            end
        end
    end
end