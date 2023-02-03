% This script was made to read data produced by smallPICOS++:
clear all
close all

% 1 - Define file names and path:
% =========================================================================
fileName = 'file_1_rank_2.h5';
dirName = '../outputFiles/HDF5_simple/';
pathName = [dirName,fileName];

fileStruct = extractStructFromHDF5(pathName);
fileStruct2 = extractStructFromHDF5('../outputFiles/HDF5/main.h5');
fileStruct3 = extractStructFromHDF5('../outputFiles/HDF5/PARTICLES_FILE_0.h5');

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
function data = extractStructFromHDF5(fileName,verbose)
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
        info=h5info(fileName,pathStr);
        
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
            data = setfield(data,fieldsName{:},varName{:},h5read(fileName,[pathStr,'/',var]));
        end
        
        % Loading attributes
        for atti=1:length(info.Attributes)
            att=info.Attributes(atti).Name;
            fields  = strsplit(pathStr,'/');   % Ignore initial blank field later
            fields(cellfun(@isempty,fields))=[];
            attName=validateFieldName(att);
            fieldsName=validateFieldName(fields);
            data = setfield(data,fieldsName{:},attName{:},h5readatt(fileName,pathStr,att));
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