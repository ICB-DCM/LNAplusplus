% runWang2010.m constructes the LNA for the model of DNA self-regulation
% introduced by Wang et al. (2010) in "Parameter inference for discretely 
% observed stochastic kinetic models using stochastic gradient descent."
% BMC Syst. Biol., 4, 99.

clear all;
close all;
clc;

% Add path
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'matlab'));

% Model:
% ======
% R1: DNA + P2  -->  DNA.P2
% R2: DNA.P2    -->  DNA + P2
% R3: DNA       -->  DNA + RNA
% R4: RNA       -->  0
% R5: P + P     -->  P2
% R6: P2        -->  P + P
% R7: RNA       -->  RNA + P
% R8: P         -->  0
    
%% Create the Matlab executable for the birth / death system
generateLNA(fullfile(fileparts(mfilename('fullpath')), 'Wang2010.xml'),'Wang2010','NONE');
% Note: For this model the symbolic calculation of the steady state fails
% as the resulting equations are already to complex.

% Add the path to the mex file
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'models', 'Wang2010'));

%% Define parameters, initial conditions and simulation time
Theta = [0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1];
species_names = {'DNA', 'DNA:P_2', 'RNA', 'P', 'P_2'};
MRE0  = [20, 1, 1, 1, 1];
Var0  = toLinear(zeros(5));
tspan = linspace(0,150,100);

%% Simulate model
[MRE,Var,sMRE,sVar,s2MRE,s2Var] = Wang2010_LNA(Theta,tspan,MRE0,Var0,0,1:5);

figure('name','Simulation');
for k = 1:5
    subplot(2,3,k);
    % Assign mean and variance for the k-th species
    m = MRE(k,:);
    s = sqrt(diag(squeeze(Var(k,k,:,:))))';
    % Plot
    h(2) = fill([tspan(1:end),tspan(end:-1:1)],...
                [m(1:end)+s(1:end),m(end:-1:1)-s(end:-1:1)],'b',...
                'EdgeColor',[0.7,0.7,1],'FaceColor',[0.7,0.7,1]); hold on;
    h(1) = plot(tspan,m,'b');
    xlabel('Time');
    ylabel('Abundance');
    title(species_names{k});
    if k == 1
        legend('mean','+/- std');
    end
end

%% Test of cross-species sensitivities
i = 5;
j = 5;
eps_theta = 1e-3;
[MRE_per,Var_per,sMRE_per,sVar_per,s2MRE_per,s2Var_per] = Wang2010_LNA(Theta+[0*[1:i-1],1,0*[i+1:8]]*eps_theta,tspan,MRE0,Var0,0,1:5);

k1 = 50;
k2 = 100;

%%
figure
subplot(1,3,1)
imagesc((Var_per(:,:,k1,k2)-Var(:,:,k1,k2))./eps_theta)
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(sVar(:,:,k1,k2,i))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc((Var_per(:,:,k1,k2)-Var(:,:,k1,k2))./eps_theta - sVar(:,:,k1,k2,i))
colorbar
title('error');

%% 2nd order test - diagonal
figure
subplot(1,3,1)
imagesc(squeeze((sVar_per(1,1,:,:,j)-sVar(1,1,:,:,j))./eps_theta))
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(squeeze(s2Var(1,1,i,j,:,:)))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc(squeeze((sVar_per(1,1,:,:,j)-sVar(1,1,:,:,j))./eps_theta) - squeeze(s2Var(1,1,i,j,:,:)))
colorbar
title('error');

%% 2nd order test - off-diagonal
figure
subplot(1,3,1)
imagesc((sVar_per(:,:,k1,k2,j)-sVar(:,:,k1,k2,j))./eps_theta)
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(s2Var(:,:,i,j,k1,k2))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc((sVar_per(:,:,k1,k2,j)-sVar(:,:,k1,k2,j))./eps_theta - s2Var(:,:,i,j,k1,k2))
colorbar
title('error');

