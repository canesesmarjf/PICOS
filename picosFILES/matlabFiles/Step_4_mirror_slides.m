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

% Extract data from HDF5 files:
% =========================================================================
% Root directory:
runType = 99;
switch runType
    case 1
        root = '/home/jfcm/Documents/compX/ARPAE/mirror_case_10/';        
    case 2
        root = '/home/jfcm/Documents/compX/ARPAE/mirror_case_11b/'; 
    case 3
        root = '/home/jfcm/Documents/compX/ARPAE/mirror_case_12/';         
    case 99
        root = '../outputFiles/';
end

% File name:
fileName = [root,'HDF5/main.h5'];

% Extract data:
extractDataFromH5;

% tt = 6; k = 1; figure; plot3(a_p{1}(:,tt),vpar_p{k}(:,tt),vper_p{k}(:,tt),'k.','markersize',3)

%% Plasma density plots:

ne = n_m{1} + n_m{2};

k = 1;
figure('color','w')
box on
hold on
hn(1) = plot(x_m,n_m{k}(:,1),'r');
hn(2) = plot(x_m,n_m{k}(:,end),'k');
plot(x_m,n_m{k}(:,1:3:end),'k')

legendText{1} = ['t = 0 [ms]'];
legendText{2} = ['t = 2 [ms]'];

hL = legend(hn,legendText);

set(hn,'lineWidth',5)
set(hL,'interpreter','latex','FontSize',13)
title('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)
ylim([0,8e19])

figure('color','w')
box on
hold on
hn(1) = plot(x_m,ne(:,1),'r');
hn(2) = plot(x_m,ne(:,end),'k');
plot(x_m,ne(:,1:3:end),'k')

legendText{1} = ['t = 0 [ms]'];
legendText{2} = ['t = 2 [ms]'];

hL = legend(hn,legendText);

set(hn,'lineWidth',5)
set(hL,'interpreter','latex','FontSize',13)
title('$n_e$ [m$^{-3}$]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)
ylim([0,8e19])

%% Computational particle density:
k = 1;
figure('color','w')
box on
hold on
hn(1) = plot(x_m,ncp_m{k}(:,1),'r');
hn(2) = plot(x_m,ncp_m{k}(:,end),'k');
% plot(x_m,ncp_m{k}(:,1:3:end),'k')

legendText{1} = ['t = 0 [ms]'];
legendText{2} = ['t = 2 [ms]'];

hL = legend(hn,legendText);

set(hn,'lineWidth',5)
set(hL,'interpreter','latex','FontSize',13)
title('$n^{cp}_m$ [m$^{-1}$]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)
% ylim([0,8e19])

%% Ion temperature:
k = 1;
figure('color','w')
subplot(1,2,1)
box on
hold on
hT(1) = plot(x_m,Tper_m{k}(:,end),'r');
hT(2) = plot(x_m,Tpar_m{k}(:,end),'k');

legendText{1} = ['t = 0 [ms]'];
legendText{2} = ['t = 2 [ms]'];

% hL = legend(hT,legendText);

set(hT,'lineWidth',5)
set(hL,'interpreter','latex','FontSize',13)
title('$T_{i\parallel}$ [eV]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)
ylim([0,3e3])

k = 2;
subplot(1,2,2)
box on
hold on
hT(1) = plot(x_m,Tper_m{k}(:,end),'r');
hT(2) = plot(x_m,Tpar_m{k}(:,end),'k');

legendText{1} = ['$T_\perp$'];
legendText{2} = ['$T_\parallel$'];

hL = legend(hT,legendText);

set(hT,'lineWidth',5)
set(hL,'interpreter','latex','FontSize',13)
title('$T_{i}$ [eV]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)
ylim([0,30e3])


%% Electric field plot:
figure('color','w')
box on
hold on
hE(1) = plot(x_m,movmean(Ex_m(:,1),10),'r');
hE(2) = plot(x_m,movmean(Ex_m(:,end),20),'k');
% plot(x_m,movmean(Ex_m(:,1:5:end),5,2),'k');

set(hE,'lineWidth',5)
title('$E_{\parallel}$ [V/m]','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)

% Electric potential:
dx = main.mesh.dx;
V = -cumtrapz(Ex_m)*dx;

figure('color','w')
box on
hold on
hV = plot(x_m,V(:,end)./Te_m,'k');
set(hV,'lineWidth',2)
ylim([0,3])

title('e$V/T_e$','Interpreter','latex','FontSize',15)
xlabel('x [m]','Interpreter','latex','FontSize',15)

%% Distribution function:
% Ion parameters:
beam.E = 25e3;
Ma = main.ions.species_1.M;

% Ion derived quantities:
Te_mean = mean(mean(Te_m));
vT = double(sqrt(2*Te_mean*e_c/Ma));

% Create mesh with ghost cells"
Nx = double(main.mesh.Nx_IN_SIM);
dx = double(main.mesh.dx);
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
for ss = 1:numel(n_m)
    fdist{ss} = fv;
end

% Calculate velocity space PDF:
% Select range:
for k = 1:numel(n_m)
rangeType = 1;
switch rangeType
    case 1
        % mirror region:
        % % ===============
        x_center = +0.0;
        x_delta  = +0.35;
        z1 = x_center - x_delta;
        z2 = x_center + x_delta;
    case 2
    case 3
    case 4  
    case 5
end

% Reference flux:
B0 = main.mesh.B0;
A0 = main.mesh.A0;
Phi0 = B0*A0;

for jj = 1:NS
    % Select the spatial range:
    rng = find(x_p{k}(:,jj) > z1 & x_p{k}(:,jj) < z2);

    Nrng = numel(rng);
        
    vvy = vper_p{k}(rng,jj);  
    vvx = vpar_p{k}(rng,jj);   
    
    [WXL,WXC,WXR,MXI,WYL,WYC,WYR,MYI] = AssignCell_2D(vvx/vT,vvy/vT,vxGrid/vT,vyGrid/vT);
    disp(['jj: ',num2str(jj)])
    
    for ii = 1:Nrng
        if ~isempty(MXI(ii))||~isempty(MYI(ii))
            ix = MXI(ii) + 2;
            iy = MYI(ii) + 2;

            vper = vper_p{k}(rng(ii),jj);
            Bm = Bx_p{k}(rng(ii),jj);
            ai = a_p{k}(rng(ii),jj);

            fv(ix - 1,iy + 0,jj) = fv(ix - 1,iy + 0,jj) + (Bm/B0)*WXL(ii)*WYC(ii)*ai./vper;
            fv(ix + 0,iy + 0,jj) = fv(ix + 0,iy + 0,jj) + (Bm/B0)*WXC(ii)*WYC(ii)*ai./vper;    
            fv(ix + 1,iy + 0,jj) = fv(ix + 1,iy + 0,jj) + (Bm/B0)*WXR(ii)*WYC(ii)*ai./vper;

            fv(ix + 0,iy - 1,jj) = fv(ix + 0,iy - 1,jj) + (Bm/B0)*WXC(ii)*WYL(ii)*ai./vper;
            fv(ix + 0,iy + 1,jj) = fv(ix + 0,iy + 1,jj) + (Bm/B0)*WXC(ii)*WYR(ii)*ai./vper;    

            fv(ix + 1,iy - 1,jj) = fv(ix + 1,iy - 1,jj) + (Bm/B0)*WXR(ii)*WYL(ii)*ai./vper;
            fv(ix + 1,iy + 1,jj) = fv(ix + 1,iy + 1,jj) + (Bm/B0)*WXR(ii)*WYR(ii)*ai./vper;

            fv(ix - 1,iy - 1,jj) = fv(ix - 1,iy - 1,jj) + (Bm/B0)*WXL(ii)*WYL(ii)*ai./vper;
            fv(ix - 1,iy + 1,jj) = fv(ix - 1,iy + 1,jj) + (Bm/B0)*WXL(ii)*WYR(ii)*ai./vper;
        end
    end

end

fdist{k} = fv;
end

%% Plot PDF:
% Mirror ratio
Bx_mirror = Bx_m(find(x_m >= 0.8,1));
% R_m =max(Bx_m(:,1))/B0;
R_m = Bx_mirror/B0;

% Loss cone angle:
theta_m = asin(1/sqrt(R_m));
theta_NBI = 50*pi/180;
vNBI = sqrt(2*e_c*25000/Ma);

% % Plotting the integral of the electric field:
figure; 
set(gcf,'Position',[1293,409,550,427])
hold on
endTime = size(Tpar_m{k},2);
rngLC = endTime-5:endTime-1;
frame=1;
ff = mean(movmean(Ex_m(:,rngLC),frame,1),2);
rng = find(x_m>x_center & x_m<z2);

delta_V = trapz(x_m(rng),ff(rng));
gg = Bx_m(:,jj);
plot(x_m(:),ff,'k');hold on;
plot(x_m(rng),ff(rng),'g');
plot(x_m(:),gg*max(ff)/max(gg),'r');
ylim([-10 10]*Te_mean);

try
%     ftotal = fdist{1};
    ftotal = fdist{1} + fdist{2};
catch
    ftotal = fdist{1};
end

% Calculating the hyperbola that defines the trapped region:
vvpar = vxGrid/vT;
vvper = sqrt((vvpar.^2 + 1*delta_V/(Te_mean))/(R_m - 1));

figure('color','w'); 
ax = gca;
rngx = 3:(NV+4-2);
rngy = 3:(NV+4-2);
Limit=12;
fmax = max(max(ftotal(:,:,end)));
for jj = 2:size(ftotal,3)-3
    contourf(ax,vxGrid/vT,vyGrid/vT,mean(ftotal(rngx,rngy,jj:jj+3),3)',10,'LineStyle','none');

    hold(ax,'on');

    lossCone(1) = line(ax,[-1,+1]*Limit,[-1,+1]*Limit*tan(theta_m));
    lossCone(2) = line(ax,[-1,+1]*Limit,[+1,-1]*Limit*tan(theta_m));
    set(lossCone,'LineWidth',2,'color','k','lineStyle','--','LineWidth',2)
        
    if 1
        NBICone(1) = line(ax,[-1,+1]*Limit,[-1,+1]*Limit*tan(theta_NBI));
        NBICone(2) = line(ax,[-1,+1]*Limit,[+1,-1]*Limit*tan(theta_NBI));
        set(NBICone,'LineWidth',2,'color','g','lineStyle','--','LineWidth',2)
        plot(ax,vNBI*cos(theta_NBI)/vT,vNBI*sin(theta_NBI)/vT,'go')
    end
    
    
    confinementRegion(1) = plot(ax,vvpar,+vvper,'k','LineWidth',1);
    confinementRegion(2) = plot(ax,vvpar,-vvper,'k','LineWidth',1);
    
    title(ax,['$f(v_\parallel,v_\perp)$ '],'interpreter','latex','FontSize',20)
    hold(ax,'off');
    axis square
    axis(ax,'image')
    xlim(ax,[-1,1]*Limit)
    ylim(ax,[0,1]*Limit)
    xlabel(ax,'$v_\parallel/v_{Te}$','interpreter','latex','FontSize',20)
    ylabel(ax,'$v_\perp/v_{Te}$','interpreter','latex','FontSize',20)
%     caxis(ax,[0,0.5]*fmax)
    colorbar(ax)
    colormap(ax,flipud(hot))
    view(ax,[0,90])
    drawnow
    pause(0.01)
end





return

%% Energy integrated over volume:
k = 1;
E_total = zeros(Nt,1);
a_total = zeros(Nt,1);
K = main.ions.species_1.K;
M = main.ions.species_1.M;

for tt = 1:Nt
    ai = a_p{k}(:,tt);
    vpar_i = vpar_p{k}(:,tt);
    vper_i = vper_p{k}(:,tt);
    v2_i = vpar_i.^2 + vper_i.^2;
    E_total(tt) = K*M*sum(ai.*v2_i)/2;
    a_total(tt) = K*sum(ai);
end

% Total energy in the system:
figure
hold on
plot(t_p*1e3,E_total)
ylim([0,2]*max(E_total))
title('total E')

% Mean energy in the system:
figure;
E_mean = E_total./a_total;
plot(t_p*1e3,E_mean/e_c)
ylim([0,2]*max(E_mean/e_c))
title('mean energy')

%% Plot ion moments:
close all

if find(strcmpi('Bx_p',vars.output) == 1) && find(strcmpi('x_p',vars.output) == 1)
    figure
    hold on
    plot(x_p{1},Bx_p{1},'k.')
    plot(x_m   ,Bx_m,'r.-')
    xlim([-2,2])
end

if find(strcmpi('n_m',vars.output) == 1)
    figure
    hold on
    for s = 1:numel(n_m)
        mesh(t_p,x_m,movmean(n_m{s},10,1));
    end
    zlim([0,1e20])
    title('n')
    ylim([-2,2])
end

if find(strcmpi('ncp_m',vars.output) == 1)
    figure
    hold on
    for s = 1:numel(ncp_m)
        mesh(t_p*1e3,x_m,movmean(ncp_m{s},10,1));
    end
    title('n_cp')
    xlabel('[ms]')
    ylim([-2,2])    
end

if find(strcmpi('Ex_m',vars.output) == 1)
    figure
    mesh(t_p,x_m,movmean(movmean(Ex_m./Te_m,10,1),4,2));
    zlim([-1,1]*30)
    caxis([-1,1]*8)
    title('Ex')
    ylim([-2,2])    
end

if find(strcmpi('u_m',vars.output) == 1)
    figure
    hold on
    for s = 1:numel(ux_m)
        Cs = sqrt( e_c*(Te_m + 3*Tpar_m{s})/main.ions.species_1.M);
        mesh(t_p,x_m,movmean(ux_m{s}./Cs,10,1));  
    end
    title('u_x/{v_T}')
    ylim([-2,2]) 
end

if find(strcmpi('Tpar_m',vars.output) == 1)
    figure
    hold on
    for s = 1:numel(Tpar_m)
        mesh(t_p,x_m,movmean(Tpar_m{s},10,1));
        Tpar_max = max(max(Tpar_m{1}));
    end
    zlim([0,5*Tpar_max])
    caxis([0,2*Tpar_max])
    title('T_par')
    ylim([-2,2])    
end

if find(strcmpi('Tper_m',vars.output) == 1)
    figure
    hold on
    for s = 1:numel(Tpar_m)    
        mesh(t_p,x_m,movmean(Tper_m{s},10,1));
    end
    zlim([0,5]*max(max(Tper_m{1})))
    title('T_per')
    ylim([-2,2])    
end

% Cross sectional area:
A0  = main.mesh.A0;
B0  = main.mesh.B0;
A_m = A0*B0./Bx_m;

if find(strcmpi('u_m',vars.output) == 1)
    % Plasma integrated flux:
    F = n_m{1}.*ux_m{1}.*A_m;

    figure
    mesh(t_p,x_m,movmean(F,10,1));
    title('particle flux')
    ylim([-2,2])    
end

%% Force balance:
% close all

% x range:
rng_x = find(x_m >= -1 & x_m <= 1);

% Plasma quantities:
M      = main.ions.species_1.M;
A0     = main.mesh.A0;
B0     = main.mesh.B0;
Te     = Te_m(rng_x,:);
Ti_par = Tpar_m{1}(rng_x,:);
Ti_per = Tper_m{1}(rng_x,:);
ne     = n_m{1}(rng_x,:);
U      = ux_m{1}(rng_x,:);
B      = Bx_m(rng_x,:);
A_m    = A0*B0./B;
x      = x_m(rng_x,:);
t      = t_p;
dx     = diff(x(1:2));
dt     = diff(t(1:2));

% Derived:
flux_density  = ne.*U;
flux          = flux_density.*A_m;
E_mean        = e_c*(Te + Ti_par + Ti_per);
power_density = flux_density.*E_mean;
power         = power_density.*A_m;
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
fr = 1;
P_par = movmean(P_par,12,1);
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
fr = 1;
F_KE_x = movmean(F_KE_x,fr,2);
F_KE_t = movmean(F_KE_t,fr,2);
F_par  = movmean(F_par ,fr,2);
F_mag  = movmean(F_mag ,fr,2);

figure; 
% rng = 26:41;
% rng = 50:60;
rng = 1:6;
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
ylim(6*[-1,1]*max(mean(F_KE_x(:,rng),2)))
xlabel('x [m]','interpreter','latex','fontSize',14)
ylabel('[Nm$^{-3}$]','interpreter','latex','fontSize',14)
grid on
try
title(['ICRF ',myDir(6:(end-7)),' [kW]'],'interpreter','latex','fontSize',14)
end

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
g = n_m{1}(:,end-3:end);
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,1e20])
title('Plasma density [m$^{-3}]$','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Computational particle density:
figure('color','w')
hold on
box on
g = ncp_m{1}(:,end-3:end);
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
% ylim([0,1e20])
title('$n_m^{cp}$ [m$^{-1}]$','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Electron tempeture:
figure('color','w')
hold on
box on
g = Te_m(:,end-3:end);
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,1.2]*fmax)
title('Electron temperature [eV]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Plasma flow:
figure('color','w')
hold on
box on
Cs = sqrt( e_c*(Te_m + 3*Tpar_m{1})/main.ions.species_1.M);
f = ux_m{1}./Cs;
Mach = mean(f(:,end-3:end),2);
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
g = n_m{1}(:,end-3:end).*Te_m(:,1)*e_c;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([0,1.2]*fmax)
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
ylim([-3,3]*1000)
title('Parallel Energy flux [MWm$^{-2}$]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

% Total power flow:
figure('color','w')
hold on
box on
g = powerFlux_par.*A_m;
f = mean(g,2);
fmax = max(f);
Bmax = max(Bx_m(:,1));
plot(x_m,movmean(f,10,1),'k','lineWidth',2);
plot(x_m,Bx_m(:,1)*fmax/Bmax)
ylim([-1,1]*6e6)
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
ylim([-1,1]*0.8*1e3)
title('Parallel convective heat flux [MWm$^{-2}$]','interpreter','latex','fontSize',14)
xlabel('x [m]','interpreter','latex','fontSize',14)

return


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
