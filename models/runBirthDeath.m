% BirthDeath.m constructes the LNA for a simple birth-death process and
% illustrates the functionality of LNA++, including the incorporation of
% measurement noise and 1st and 2nd order senstivity analysis.

clear all;
close all;
clc;

% Add path
addpath('../matlab')

% Model of bearth-death process:
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
generateLNA('BirthDeath/BirthDeath.xml', 'BirthDeath', 'BOTH');

% add the path to the mex file
addpath('BirthDeath/');

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
title('Macroscopic mean of protein');
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
ObsIndex = [1,2]; % observable is second species (= mRNA,protein)
VarNoise = [10,50]; % variance of measurement noise

% solve LNA and compute measured distribution
[MRE,Var] = BirthDeath_LNA(Theta,tspan,ObsIndex,VarNoise);

% plot results
figure('Name','Simulation: IC = steady state; observable = mRNA and protein')
% - mean of mRNA
subplot(2,2,1);
plot(tspan,MRE(1,:));
xlabel('Time');
ylabel('Protein');
title('Macroscopic mean of protein');
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
[MRE,Var] = BirthDeath_LNA(Theta,tspan,ObsIndex,VarNoise);

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

MRE0 = [2,200]; % inital values (E(mRNA),E(Protein))
Var0 = toLinear([0,0;0,0]); % initia co-variances (cov[mRNA,mRNA],cov(mRNA,Protein),cov(Protein,Protein))

% solve LNA and compute measured distribution
[MRE,Var] = BirthDeath_LNA(Theta,tspan,ObsIndex,VarNoise,MRE0,Var0);

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
[MRE,Var,Sens_MRE,Sens_Var] = BirthDeath_LNA(Theta,tspan,ObsIndex,VarNoise,MRE0,Var0);

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
[MRE,Var,Sens_MRE,Sens_Var,Sens2_MRE,Sens2_Var] = BirthDeath_LNA(Theta,tspan,ObsIndex,VarNoise,MRE0,Var0);

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