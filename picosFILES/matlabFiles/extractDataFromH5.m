disp('Begin data extraction ...')
disp('********************************************************************')
disp('')

% Read MAIN file:
% =========================================================================

% Get data:
main = H5_to_struct(fileName);

% Number of ranks for fields:
ranksFields = double(main.numMPIsFields);

% Number of ranks for particles:
ranksParticles = double(main.numMPIsParticles);

% Total number of mesh points:
Nx = main.mesh.Nx_IN_SIM;

% Mesh points per MPI:
NxPerMPI = double(main.mesh.Nx_PER_MPI);

% x axis:
x_m  = double(main.mesh.x_m);

% Number of ions species:
numIonSpecies = main.ions.numberOfParticleSpecies;

% Ion species parameters:
for ii = 1:numIonSpecies
    fieldName = ['species_',num2str(ii)];
    ionParameters{ii} = main.ions.(fieldName);
end

% Get number of time outputs:
% =========================================================================
fileName  = [root,'HDF5/FIELDS_FILE_0.h5'];
fields    = H5_to_struct(fileName);
timeSteps = fieldnames(fields);
Nt = length(timeSteps);
timeMax = double(fields.(['t',num2str(Nt-1)]).time);


% Determine the outputted variables:
% =========================================================================
% fields:
vars.fields = fieldnames(fields.t0.fields);

% plasma:
fileName  = [root,'HDF5/PARTICLES_FILE_0.h5'];
plasma      = H5_to_struct(fileName);
vars.ions      = fieldnames(plasma.t0.ions.species_1);
vars.electrons = fieldnames(plasma.t0.electrons);

vars.output = [vars.fields;vars.ions;vars.electrons];

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
    fileName = [root,'HDF5/FIELDS_FILE_',num2str(rr-1),'.h5'];
    
    % Get data:
    fields = H5_to_struct(fileName);
    
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
    N_CP_MPI = double(ionParameters{1}.N_CP_MPI_Output);

    % Total number of particles:
    N_CP = N_CP_MPI*ranksParticles;
        
    % Initialize data:
    t_p        = zeros(Nt ,1 );
    
    % Particle states:
    if sum(strcmpi('x_p',vars.output))
        x_p{ss}    = zeros(N_CP,Nt);
    end
    
    if sum(strcmpi('v_p',vars.output))
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
        Te_m = zeros(Nx,Nt);
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
    fileName = [root,'HDF5/PARTICLES_FILE_',num2str(rr-1),'.h5'];
    
    % Get data:
    plasma = H5_to_struct(fileName);
    
    % Loop over all ion species:
    for ss = 1:numIonSpecies
        
        % Species field:
        speciesField = ['species_',num2str(ss)];
              
        % Number of particles per rank:
        N_CP_MPI = double(ionParameters{ss}.N_CP_MPI_Output);
                        
        % Data range:
        rng_p = 1 + N_CP_MPI*(rr-1) : 1 : N_CP_MPI*rr;

        % Loop over all time:
        for tt = 1:Nt
            % Obtain time stamp:
            kk = str2double(timeSteps{tt}(2:end)) + 1;
            
            % Time:
            t_p(kk) = plasma.(timeSteps{tt}).time;
            
           % Particle states:
            if sum(strcmpi('x_p',vars.output))
                x_p{ss}(rng_p,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).x_p;
            end
            
            if sum(strcmpi('v_p',vars.output))
                vpar_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).v_p(:,1);
                vper_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).v_p(:,2);   
            end
            
            if sum(strcmpi('a_p',vars.output))
                a_p{ss}(rng_p,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).a_p;
            end
            
            if sum(strcmpi('mu_p',vars.output))
                mu_p{ss}(rng_p,kk)   = plasma.(timeSteps{tt}).ions.(speciesField).mu_p;
            end        
        
            % Particle-defined moments:
            if sum(strcmpi('n_p',vars.output))
                n_p{ss}(rng_p,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).n_p;
            end
            
            if sum(strcmpi('nv_p',vars.output))
                nv_p{ss}(rng_p,kk)   = plasma.(timeSteps{tt}).ions.(speciesField).nv_p;
            end
            
            if sum(strcmpi('Tpar_p',vars.output))
                Tpar_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).Tpar_p;
            end
            
            if sum(strcmpi('Tper_p',vars.output))
                Tper_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).Tper_p;
            end
            
            if sum(strcmpi('Te_p',vars.output))
                Te_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).Te_p;
            end
            
            if sum(strcmpi('u_p',vars.output))
                ux_p{ss}(rng_p,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).u_p.x;               
            end
    
            % Mesh-defined moments:
            if rr == 1
                if sum(strcmpi('n_m',vars.output))
                    n_m{ss}(:,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).n_m;
                end

                if sum(strcmpi('ncp_m',vars.output))
                    ncp_m{ss}(:,kk)    = plasma.(timeSteps{tt}).ions.(speciesField).ncp_m;
                end
                
                if sum(strcmpi('nv_m',vars.output))
                    nv_m{ss}   = plasma.(timeSteps{tt}).ions.(speciesField).nv_m;
                end
                
                if sum(strcmpi('Tpar_m',vars.output))
                    Tpar_m{ss}(:,kk) = plasma.(timeSteps{tt}).ions.(speciesField).Tpar_m;
                end
                
                if sum(strcmpi('Tper_m',vars.output))
                    Tper_m{ss}(:,kk) = plasma.(timeSteps{tt}).ions.(speciesField).Tper_m;
                end
                
                if sum(strcmpi('Te_m',vars.output))
                    Te_m(:,kk) = plasma.(timeSteps{tt}).electrons.Te_m;
                end
                
                if sum(strcmpi('u_m',vars.output))
                    ux_m{ss}(:,kk)   = plasma.(timeSteps{tt}).ions.(speciesField).u_m.x;                
                end
            end
                
            % Particle-defined fields:
            if sum(strcmpi('Bx_p',vars.output))
                Bx_p{ss}(rng_p,kk)   = plasma.(timeSteps{tt}).ions.(speciesField).Bx_p;
            end
            
            if sum(strcmpi('dBx_p',vars.output))
                dBx_p{ss}(rng_p,kk)  = plasma.(timeSteps{tt}).ions.(speciesField).dBx_p;
            end
            
            if sum(strcmpi('ddBx_p',vars.output))
                ddBx_p{ss}(rng_p,kk) = plasma.(timeSteps{tt}).ions.(speciesField).ddBx_p;
            end
            
            if sum(strcmpi('Ex_p',vars.output))
                Ex_p{ss}(rng_p,kk)   = plasma.(timeSteps{tt}).ions.(speciesField).Ex_p;
            end         
                                    
        end
        
    end
end
% cd ..

disp('Particle data extraction completed!')
disp('********************************************************************')
disp('')