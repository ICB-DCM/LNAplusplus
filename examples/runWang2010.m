% runWang2010.m constructes the LNA for the model of DNA self-regulation
% introduced by Wang et al. 2010.

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

% add the path to the mex file
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'models', 'Wang2010'));

%% Define parameters, initial conditions and simulation time
Theta = zeros(1,8);%[0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1];
%Theta = [0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1];
species_names = {'DNA', 'DNAP2', 'RNA', 'P', 'P2'};
MRE0  = [20, 0, 0, 0, 0];
Var0  = toLinear(zeros(5));
tspan = linspace(0,150,100);

%% Simulate model
[MRE,Var,sMREsdVar] = Wang2010_LNA(Theta,tspan,MRE0,Var0,0,1:5);

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