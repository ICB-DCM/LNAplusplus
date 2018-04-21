% LinearChain.m constructes the LNA for a linear reaction chain of varying
% length. This is used to assess the scalability of the symbolic
% construction as well as the numerical simulation without sensitivities,
% with 1st order sensitivities and with 2nd order sensitivities.

clear all;
close all;
clc;

% Add path
addpath('../matlab')
addpath('./linearChain')

% Model of bearth-death process:
% ==============================
% X_{i} -> X_{i+1} for i = 1 to N-1

%% Define evaluation parameters
R = 10; % number of simulation runs
N_max = 6; % maximum chain length (the minimum is 3)

%% Initialize matrices for storing computation and run time
computeTimeMat = nan(N_max-2,1); % time for model construction and compilation
runTimeMat0    = nan(N_max-2,R); % time for simulation w/o sensitivities
runTimeMat1    = nan(N_max-2,R); % time for simulation with 1st order sensitivities
runTimeMat2    = nan(N_max-2,R); % time for simulation with 2nd order sensitivities

%% estimate compute and run time for models
for N = 3:N_max
    syms t; % time
    phi   = sym('phi'  ,[1,  N]); % symbolic variables for macroscopic mean
    Theta = sym('Theta',[1,2*N]); % symbolic variables for paraeters
    assume(Theta, 'positive')
    
    % Generate stochiometic matric and flux vector
    SS = S_gen_linearChain(N);
    F = @(phi,t,Theta) [Theta(1),Theta(2)*phi(1),phi(1:end-1).*Theta(3:length(phi)+1),phi(2:end).*Theta(length(phi)+2:end)].';
    FF = matlabFunction(F(phi,t,Theta),'vars',{phi, t, Theta});

    % Generate model
    tic
    modelName = ['chain' int2str(N)]; % name
    generateLNAComponents(modelName, SS, FF, phi, Theta, 'NONE')
    compileLNA(modelName); 
    computeTimeMat(N-2) = toc;
    
    % Define function handle
    fHandle = eval(['@chain' int2str(N) '_LNA']);
    
    % Add path
    addpath(modelName);
    
    % Run simulation
    disp('Running simulations...')
    tspan = 0:0.1:2.5;
    for i = 1:R
        % Generate parameters and initial conditions
        Theta = rand(1,2*N)*100;
        MRE0  = rand(1,N)*1000;
        Var0  = zeros(1,N*(N+1)/2);
        
        % Simulate w/o sensitivities
        tic
        [M,S] = fHandle(Theta,tspan,1:N,0,MRE0,Var0);
        runTimeMat0(N-2,i) = toc;        
        fprintf('.')

        % Simulate with 1st order sensitivities
        tic
        [M1,S1,dM,dS] = fHandle(Theta,tspan,1:N,0,MRE0,Var0);
        runTimeMat1(N-2,i) = toc;        
        fprintf('.')
        
        % Simulate with 2nd order sensitivities
        tic
        [M2,S2,dM,dS,d2M,d2S] = fHandle(Theta,tspan,1:N,0,MRE0,Var0);
        runTimeMat2(N-2,i) = toc;        
        fprintf('.')
    end
    fprintf('\n')
end

%% Plot results
figure('Name','Computation time for model construction');
semilogy(3:N_max, computeTimeMat,'.-','Color',[0.4 0.4 0.4],'MarkerSize',30);
xlabel('Network size');
ylabel('Model construction time (s)');

figure('Name','Computation time for model simulation');
semilogy(3:N_max,median(runTimeMat0'),'k.-','MarkerSize',30); hold on;
semilogy(3:N_max,median(runTimeMat1'),'b.-','MarkerSize',30);
semilogy(3:N_max,median(runTimeMat2'),'r.-','MarkerSize',35);
xlabel('Network size');
ylabel('Median run time (s)');
legend('w/o sensitivities','1st order sensitivities','2nd order sensitivities');
