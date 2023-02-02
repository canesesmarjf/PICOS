% This script was made to read data produced by smallPICOS++:
clear all
close all

% 1 - Define file names and path:
% =========================================================================
fileName = 'file_1_rank_2.h5';
dirName = '../outputFiles/HDF5_simple/';
pathName = [dirName,fileName];

% 2- Get group names:
% =========================================================================
info = h5info(pathName);
for i = 1:numel(info.Groups)
    groupNames{i} = info.Groups(i).Name;
    for j = 1:numel(info.Groups(i).Datasets)
        datasetNames{i}{j} = info.Groups(i).Datasets(j).Name;
        info.Groups(i).Datasets(j).Name
    end
end

% 3 - Read data:
% =========================================================================
for i = 1:numel(info.Groups)
    group = groupNames{i}(2:end);
    disp(["Group = ",group]);
    for j = 1:numel(info.Groups(i).Datasets)
        var = info.Groups(i).Datasets(j).Name;
        values = h5read(pathName,[groupNames{i},'/',datasetNames{i}{j}]);
        
        d.(group).(var) = values;
        data{i}{j} = values;
    end
end 

% Assign data:
% =========================================================================
for i = 1:numel(info.Groups)
    a_p{i} = data{i}{1};
    v_p{i} = data{i}{2};
    x_p{i} = data{i}{3};
end 

% Calculate histograms for each species:
% =========================================================================
x_m = d.fields.x_m;
dx = x_m(2) - x_m(1);
x_edges = x_m(1:end-1) + dx/2;
[counts_0,~] = histcounts(d.ions_0.x_p,x_edges,"Normalization","count");
[counts_1,~] = histcounts(d.ions_1.x_p,x_edges,"Normalization","count");

% Plot data:
% =========================================================================
figure('color','w')
subplot(2,1,1)
box on
hold on
plot(d.fields.x_m,d.ions_0.ncp_m*dx,'k.-')
plot(d.fields.x_m,d.ions_1.ncp_m*dx,'r.-')

plot(x_m(2:end-1),counts_0,'k','LineWidth',2)
plot(x_m(2:end-1),counts_1,'r','LineWidth',2)

ylim([0,1.2]*max(d.ions_0.ncp_m*dx))

subplot(2,1,2)
box on
hold on
plot(d.fields.x_m,d.ions_0.n_m,'k.-')
plot(d.fields.x_m,d.ions_1.n_m,'r.-')
ylim([0,1.2]*max(d.ions_0.n_m))
xlabel(['$x$ [m]'],'interpreter','latex','FontSize',18)
title(['Real particle density $n_m^{R}$ [m$^{-3}$]'],'interpreter','latex','FontSize',18)

figure('color','w')
box on
hold on
plot(d.fields.x_m,d.ions_0.Tpar_m,'k.-')
plot(d.fields.x_m,d.ions_1.Tpar_m,'r.-')
plot(d.fields.x_m,d.ions_0.Tper_m,'k.-')
plot(d.fields.x_m,d.ions_1.Tper_m,'r.-')
ylim([0,+1]*1.2*max(d.ions_0.Tper_m))
xlabel(['$x$ [m]'],'interpreter','latex','FontSize',18)
title(['Ion temperature $T_i$ [eV]'],'interpreter','latex','FontSize',18)

figure('color','w')
box on
hold on
plot(d.fields.x_m,d.ions_0.upar_m,'k.-')
plot(d.fields.x_m,d.ions_1.upar_m,'r.-')
ylim([-1,+1]*1.2*max(d.ions_0.upar_m))
xlabel(['$x$ [m]'],'interpreter','latex','FontSize',18)
title(['Ion temperature $T_i$ [eV]'],'interpreter','latex','FontSize',18)

figure
plot(d.ions_0.x_p,d.ions_0.v_p(:,1),'k.');

Npop = numel(d.ions_0.a_p);
Nsample = round(Npop*0.6);
R1 = randperm(Npop);
rng = R1(1:Nsample);

[f_pop,~] = histcounts2(d.ions_0.x_p     ,d.ions_0.v_p(:,1)  ,64,'Normalization','pdf');
[f_sam,~] = histcounts2(d.ions_0.x_p(rng),d.ions_0.v_p(rng,1),64,'Normalization','pdf');

figure; 
surf(f_pop,'LineStyle','none')
figure; 
surf(f_sam,'LineStyle','none')

figure;
subplot(2,1,1)
hold on
plot(trapz(f_pop,1))
plot(trapz(f_sam,1))

subplot(2,1,2)
hold on
plot(trapz(f_pop,2))
plot(trapz(f_sam,2))

%% Functions: