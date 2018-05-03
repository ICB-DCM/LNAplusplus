% runBirthDeath.m constructs the LNA for a simple birth-death process and
% illustrates the functionality of LNA++, including the incorporation of
% measurement noise and 1st and 2nd order senstivity analysis.

clear all;
close all;
clc;

% Add path
exampleDir = fullfile(fileparts(mfilename('fullpath')));
lnaDir = fullfile(fileparts(fileparts(mfilename('fullpath'))));
addpath([lnaDir '/matlab']);

% Model of birth-death process:
% ==============================
%          k_m
% DNA      -->  DNA+mRNA
%          g_m
% mRNA     -->  0
%          k_p
% mRNA     -->  mRNA + protein
%          g_p
% protein  -->  0

%% Create the Matlab executable for the birth / death system
generateLNA([exampleDir '/BirthDeath.xml'], 'BirthDeath', 'BOTH');

% add the path to the mex file
addpath([lnaDir '/models/BirthDeath/']);

%% Setting
Theta = [20,25,10,1]; % parameters
Theta_name = {'k_m','k_p','g_m','g_p'}; % parameter names (for visualization)
tspan = 0:0.1:10; % observation times

%% Simulate: IC = steady state
% solve LNA
[MRE,Var] = BirthDeath_LNA(Theta,tspan);

% plot results
figure('Name','Simulation: IC = Steady state')
% - mean of mRNA
subplot(2,2,1);
plot(tspan,MRE(1,:));
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of mRNA');
% - mean of protein
subplot(2,2,3);
plot(tspan,MRE(2,:));
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of protein');
% - autocovariance of mRNA
subplot(2,2,2);
imagesc(tspan,tspan,squeeze(Var(1,1,:,:)));
xlabel('Time');
ylabel('Time');
title('Autocovariance of mRNA');
% - autocovariance of Protein
subplot(2,2,4);
imagesc(tspan,tspan,squeeze(Var(2,2,:,:)));
xlabel('Time');
ylabel('Time');
title('Autocovariance of protein');

%% Simulate: IC = steady state; observable = mRNA and protein
ObsIndex = [1,2]; % Both species are observed (= mRNA,protein)
VarNoise = [10,50]; % variance of measurement noise

% solve LNA and compute measured distribution
[MRE,Var] = BirthDeath_LNA(Theta,tspan,[],[],VarNoise,ObsIndex);
% plot results
figure('Name','Simulation: IC = steady state; observable = mRNA and protein')
% - mean of mRNA
subplot(2,2,1);
plot(tspan,MRE(1,:));
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of mRNA');
% - mean of protein
subplot(2,2,3);
plot(tspan,MRE(2,:));
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of protein');
% - autocovariance of mRNA
subplot(2,2,2);
imagesc(tspan,tspan,squeeze(Var(1,1,:,:)));
xlabel('Time');
ylabel('Time');
title('Autocovariance of mRNA');
% - autocovariance of Protein
subplot(2,2,4);
imagesc(tspan,tspan,squeeze(Var(2,2,:,:)));
xlabel('Time');
ylabel('Time');
title('Autocovariance of protein');

%% Simulate: IC = steady state; observable = protein
ObsIndex = 2; % observable is second species (= protein)
VarNoise = 50; % variance of measurement noise

% solve LNA and compute measured distribution
[MRE,Var] = BirthDeath_LNA(Theta,tspan,[],[],VarNoise,ObsIndex);

% plot results
figure('Name','Simulation: IC = steady state; observable = protein')
% - mean of protein
subplot(1,2,1);
plot(tspan,MRE);
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of protein');
% - autocovariance of protein
subplot(1,2,2);
imagesc(tspan,tspan,Var);
xlabel('Time');
ylabel('Time');
title('Autocovariance of protein');

%% Simulate: IC = no steady state; observable = protein
ObsIndex = 2; % observable is second species (= protein)
VarNoise = 0; % variance of measurement noise

MRE0 = [2,200]; % initial values (E(mRNA),E(Protein))
Var0 = toLinear([0,0;0,0]); % initial co-variances (cov[mRNA,mRNA], cov(mRNA,Protein), cov(Protein,Protein))

% solve LNA and compute measured distribution
[MRE,Var] = BirthDeath_LNA(Theta,tspan,MRE0,Var0,VarNoise,ObsIndex);

% plot results
figure('Name','Simulation: IC = no steady state; observable = protein')
% - mean of protein
subplot(1,2,1);
plot(tspan,MRE);
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of protein');
% - autocovariance of protein
subplot(1,2,2);
imagesc(tspan,tspan,Var);
xlabel('Time');
ylabel('Time');
title('Autocovariance of protein');

%% Simulate: IC = no steady state; observable = protein; sensitivity = 1st order
% solve LNA and compute measured distribution
[MRE,Var,Sens_MRE,Sens_Var] = BirthDeath_LNA(Theta,tspan,MRE0,Var0,VarNoise,ObsIndex);

% plot results
figure('Name','Simulation: IC = no steady state; observable = protein; sensitivity = 1st order')
for i=1:length(Theta)
    % - mean of protein
    subplot(2,length(Theta),i);
    plot(tspan,Sens_MRE(i,:)); hold on;
    ylabel('Sensitivity of mean')
    xlabel('Time')
    title(Theta_name{i});
    % - autocovariance of protein
    subplot(2,length(Theta),length(Theta)+i);
    imagesc(tspan,tspan,Sens_Var(:,:,i));
    xlabel('Time');
    ylabel('Time');
    title(Theta_name{i});
    h = colorbar;
    h.Label.String = 'Sensitivity of autocovariance';
end

%% Simulate: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order
% solve LNA and compute measured distribution
[MRE,Var,Sens_MRE,Sens_Var,Sens2_MRE,Sens2_Var] = BirthDeath_LNA(Theta,tspan,MRE0,Var0,VarNoise,ObsIndex);

% plot results
figure('Name','Simulation of mean: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order')
for i=1:length(Theta)
    for j=1:length(Theta)
        % - mean of protein
        subplot(length(Theta),length(Theta),length(Theta)*(i-1)+j);
        plot(tspan, squeeze(Sens2_MRE(i,j,:)), 'LineWidth', 3)
        ylabel('Sensitivity of mean')
        xlabel('Time')
        title(['(' Theta_name{i} ',' Theta_name{j} ')']);
    end
end

figure('Name','Simulation of autocovariance: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order')
for i=1:length(Theta)
    for j=1:length(Theta)
        % - autocovariance of protein
        subplot(length(Theta),length(Theta),length(Theta)*(i-1)+j);
        imagesc(tspan,tspan,squeeze(Sens2_Var(i,j,:,:)));
        xlabel('Time');
        ylabel('Time');
        title(['(' Theta_name{i} ',' Theta_name{j} ')']);
        h = colorbar;
        h.Label.String = 'Sensitivity';
    end
end


%% Test of cross-species sensitivities
i = 4;
j = 4;
eps_theta = 1e-4;
[MRE,Var,sMRE,sVar,s2MRE,s2Var] = BirthDeath_LNA(Theta,tspan);
[MRE_per,Var_per,sMRE_per,sVar_per,s2MRE_per,s2Var_per] = BirthDeath_LNA(Theta+[0*[1:i-1],1,0*[i+1:4]]*eps_theta,tspan);

k1 = 50;
k2 = 100;

% 1st order sensitivity matrix
figure('Name','Test of 1st order sensitivity matrix')
subplot(1,3,1)
imagesc((Var_per(:,:,k1,k2)-Var(:,:,k1,k2))./eps_theta)
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(sVar(:,:,i,k1,k2))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc((Var_per(:,:,k1,k2)-Var(:,:,k1,k2))./eps_theta - sVar(:,:,i,k1,k2))
colorbar
title('error');

% 2nd order sensitivity matrix for temporal cross-covariance of protein (species 2) abundance
figure('Name','Test of 2nd order sensitivity matrix for temporal cross-covariance of protein abundance')
subplot(1,3,1)
imagesc(squeeze((sVar_per(2,2,j,:,:)-sVar(2,2,j,:,:))./eps_theta))
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(squeeze(s2Var(2,2,i,j,:,:)))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc(squeeze((sVar_per(2,2,j,:,:)-sVar(2,2,j,:,:))./eps_theta) - squeeze(s2Var(2,2,i,j,:,:)))
colorbar
title('error');

% 2nd order sensitivity matrix for temporal cross-covariance of two time points
figure('Name','Test of 2nd order sensitivity matrix for temporal cross-covariance of two time points')
subplot(1,3,1)
imagesc((sVar_per(:,:,j,k1,k2)-sVar(:,:,j,k1,k2))./eps_theta)
colorbar
title('finite differences');
subplot(1,3,2)
imagesc(s2Var(:,:,i,j,k1,k2))
colorbar
title('analytical sensitivities');
subplot(1,3,3)
imagesc((sVar_per(:,:,j,k1,k2)-sVar(:,:,j,k1,k2))./eps_theta - s2Var(:,:,i,j,k1,k2))
colorbar
title('error');

