% Simple birth/death model
% (C) Justin Feigelman
% justin.feigelman@helmholtz-muenchen.de
% Helmholtz Zentrum Muenchen
% 2014

%%%%%%%%%%%%%%% create the Matlab executable for the birth / death system %%%%%%%%%%%%%%% 
%% define the stoichiometric matrix
clear all

    model = 'BirthDeath';

% define the stoichiometric matrix
% the model correspond to the system
%     k_m
% DNA -->  DNA+mRNA
%      g_m
% mRNA --> 0
%       k_p
% mRNA  -->  mRNA + protein
%         g_p
% protein -->  0

S = [1 -1 0 0; % change to mRNA
    0 0 1 -1]; % change to protein

addpath('../matlab')
%  define symbolic variables
syms k_m k_p g_m g_p real
syms m p real

% define reaction fluxes
f = @(phi,t,Theta) ...
    [ Theta(1), Theta(3)*phi(1), Theta(2)*phi(1), Theta(4)*phi(2)];

phi     = [m,p]; % state vector
Theta   = [k_m, k_p, g_m, g_p]; % parameter vector
npar 	= length(Theta); % 4 parameters

% generate the C code for the components of the LNA model
generateLNAComponents(model, S, f, phi, Theta, 'BOTH')

% compile code
compileLNA(model, S, npar); 

%%%%%%%%%%%%%%% run some simulations and show the results %%%%%%%%%%%%%%%     
%% set up the simulations

% define model parameters, in same order as above
Theta = [20, 25, 10, 1];

% define observation times for the simulation
tspan = 0:0.1:10;

% add the path to the mex file
addpath('BirthDeath/')

%% compute solution assuming steady state and no measurement error

% compute the solution
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan);
toc

% plot MRE
close all
figure('Name','BirthDeath LNA Simulation: Steady State')
subplot(121)
plot(tspan,MRE(2,:),'LineWidth',2) % just plot protein
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')

% plot variance
subplot(122)
imagesc(squeeze(Var(2,2,:,:)))
xlabel('Time')
ylabel('Time')
title('Protein Autocovariance matrix')

%% specify an observed variable (Protein)
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2); % outputs protein (species 2) only
toc

%% specify measurement error (same for each species)
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2, 50); % mRNA and protein have measurement variance 50
toc

figure('Name', 'BirthDeath LNA Simulation: Added Measurement Error')
subplot(121)
plot(tspan,MRE,'LineWidth',2)
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, Var)
xlabel('Time')
ylabel('Time')
title('Protein Autocovariance matrix')

%% specify measurement error (each separately)
close all
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, [1,2], [50, 100]); % measurement error 50 and 100 for mRNA and protein, respectively
toc

figure('Name','BirthDeath LNA Simulation: Added Measurement Error for All Species')
subplot(121)
plot(tspan,MRE(2,:),'LineWidth',2)
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, squeeze(Var(2,2,:,:)))
xlabel('Time')
ylabel('Time')
title('Protein Autocovariance matrix')

%% specify initial values, variances

y0 = [2, 200]; % inital values [mRNA, Protein]
v0 = [0, 0, 0 ]; % variances

tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
toc

figure('Name', 'BirthDeath LNA Simulation: Specify Initial Conditions')
subplot(121)
plot(tspan,MRE,'LineWidth',2)
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, Var)
xlabel('Time')
ylabel('Time')
title('Protein Autocovariance matrix')

%% Compute 1st order sensitivities
tic
[MRE, Var, Sens_MRE, Sens_Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
toc
%% Plot first order sensitivity plots: analytical
hSens1_MRE = figure('Name','BirthDeath First Order Sensitivities of the MRE');
%% first order sensitivity plots: analytical
close all
figure(1)
set(gcf, 'PaperPosition', [0.25, 2.5, 6,6])
labels = {'k_m','k_p','g_m','g_p'};
for i=1:4
    subplot(2,2,i)
    plot(tspan, Sens_MRE(i,:), 'LineWidth', 2)
    title(['Sensitivity ' labels{i}])
    xlabel('Time')
    set(gca,'FontSize',20)    
end

hSens1_Var = figure('Name','BirthDeath First Order Sensitivities of the Autocovariance');
%% First order sensitivity of the variance
figure(2)
set(gcf, 'PaperPosition', [0.25, 2.5, 8,6])
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan,Sens_Var(:,:,i))
    title(['Sensitivity ' labels{i}])
%     xlabel('Time (a.u.)')
%     ylabel('Time (a.u.)')
    set(gca,'FontSize',20)
    colorbar
end
% colorbar
saveas(gcf, '~/Projects/LNA++/paper/FirstOrderVarSensitivityAnalytical.pdf')


%% compute finite difference approximation to first order sensitivities
delta=0.01;
clear FD_MRE FD_Var

for i=1:4
    [MREa, Vara] = BirthDeath_LNA(Theta, tspan, 2,0, y0, v0);
    Thetab = Theta;
    Thetab(i) = Thetab(i)*(1+delta);
    [MREb, Varb] = BirthDeath_LNA(Thetab, tspan, 2,0, y0, v0);    
    
    FD_MRE(i,:)     = (MREb-MREa)/(Theta(i)*delta);
    FD_Var(i,:,:)   = (Varb-Vara)/(Theta(i)*delta);
end


%% comparison of the MRE first order sensitivities
figure(hSens1_MRE)
for i=1:4
    subplot(2,2,i)
    hold all
    plot(tspan, FD_MRE(i,:), 'o', 'MarkerSize', 2)
%     title(['Sensitivity ' labels{i}])
    xlabel('Time (a.u.)')
    set(gca,'FontSize',20)
end
legend('NorthEast', {'LNA++','F.D.'})
%%
saveas(gcf, '~/Projects/LNA++/paper/FirstOrderMRESensitivityComparison.pdf')
% print('~/Projects/LNA++/paper/FirstOrderMRESensitivityComparison.eps','-depsc')
saveas(gcf, '~/Projects/LNA++/paper/FirstOrderMRESensitivityComparison.fig')

%% plot finite difference first order var sensitivity

figure('Name', 'BirthDeath F.D. First Order Sensitivities of the Autocovariance')
set(gcf, 'PaperPosition', [0.25, 2.5, 8,6])
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan,squeeze(FD_Var(i,:,:)))
    title(['Sensitivity ' labels{i}])
%     xlabel('Time (a.u.)')
%     ylabel('Time (a.u.)')
    set(gca,'FontSize',20)
    colorbar
end

saveas(gcf, '~/Projects/LNA++/paper/FirstOrderVarSensitivityFD.pdf')

%% compute second order sensitivities

tic
[MRE, Var, Sens_MRE, Sens_Var, Sens2_MRE, Sens2_Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
toc
%% compute second order sensitivities assuming steady state
tic
[MRE, Var, Sens_MRE, Sens_Var, Sens2_MRE_SS, Sens2_Var_SS] = BirthDeath_LNA(Theta, tspan, 2, 0);
toc
%% plot second order sensitivity of the MRE

hSens2_MRE = figure('Name', 'BirthDeath Second Order Sensitivities of the MRE');
labels = {'k_m','k_p','g_m','g_p','m_0','p_0'};
set(gcf,'PaperPosition',[0.25,2.5,6,6])
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        plot(tspan, squeeze(Sens2_MRE(i,j,:)), 'LineWidth', 2)

%         xlabel('Time (a.u.)')
        if i==1
            title(labels{j})
        end        
        if j==1
            ylabel(labels{i}, 'FontWeight', 'bold');
        end
        set(gca,'FontSize',14)

        k=k+1;
    end
end
%% plot second order sensitivity of the Autocovariance
close all
figure('Name', 'BirthDeath Second Order Sensitivities of the Autocovariance')
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
                
        clims(i,j,1) = min(min(squeeze(Sens2_Var(i,j,:,:))));
        clims(i,j,2) = max(max(squeeze(Sens2_Var(i,j,:,:))));
        
        imagesc(tspan,tspan,squeeze(Sens2_Var(i,j,:,:))) % squeeze(clims(i,j,:))')
%         title(sprintf('Sensitivity(%s,%s)', labels{i}, labels{j}))
        if i==1
            title(labels{j})
        end        
%         if j==1
%             ylabel({labels{i},'Time (a.u.)'}, 'FontWeight', 'bold');
%         end
        set(gca,'FontSize',14)

%         ylabel('Time (a.u.)')
%         xlabel('Time (a.u.)','FontWeight','bold')
        k=k+1;
        colorbar
    end
end

saveas(gcf,'~/Projects/LNA++/paper/Sens2Var.pdf')


%% compute second order finite difference 
clear FD_MRE FD_Var
% Theta = [10 10 2 1];
% y0 = [5 50];
Theta = [20 25 10 1];
y0 = [2 200];
% v0 = [10,1,500];
merr=0;
[MREa, Vara] = BirthDeath_LNA(Theta, tspan, 2, merr, y0, v0);
 
% deltaVec = 10.^-[1:5]; % produces reasonable results with this perturbation size
deltaVec=10^-1.3;
for d=1:length(deltaVec)
    delta=deltaVec(d);

    tic
    for i=1:4
        for j=1:4
            fprintf('.')
            if i==j
                Theta_f = Theta; Theta_b=Theta;
                Theta_f(i) = Theta(i)*(1+delta);
                Theta_b(i) = Theta(i)*(1-delta);

                [MRE_f, Var_f] = BirthDeath_LNA(Theta_f, tspan, 2, merr, y0, v0);
                [MRE_b, Var_b] = BirthDeath_LNA(Theta_b, tspan, 2, merr, y0, v0);    

                FD_MRE(i,i,d,:)     = (MRE_f+MRE_b-2*MREa)/(Theta(i)*delta)^2;
                FD_Var(i,i,d,:,:)   = (Var_f+Var_b-2*Vara)/(Theta(i)*delta)^2;
            else
                Theta_ff = Theta; Theta_f0 = Theta; Theta_0f = Theta; 
                Theta_b0 = Theta; Theta_0b = Theta; Theta_bb = Theta;

                Theta_ff([i j]) = Theta([i j])*(1+delta);           
                Theta_f0(i) = Theta(i)*(1+delta);
                Theta_0f(j) = Theta(j)*(1+delta);

                Theta_b0(i) = Theta(i)*(1-delta);
                Theta_0b(j) = Theta(j)*(1-delta);
                Theta_bb([i j]) = Theta([i j])*(1-delta);

                [MRE_ff, Var_ff] = BirthDeath_LNA(Theta_ff, tspan, 2, merr, y0, v0);
                [MRE_f0, Var_f0] = BirthDeath_LNA(Theta_f0, tspan, 2, merr, y0, v0);
                [MRE_0f, Var_0f] = BirthDeath_LNA(Theta_0f, tspan, 2, merr, y0, v0);                

                [MRE_bb, Var_bb] = BirthDeath_LNA(Theta_bb, tspan, 2, merr, y0, v0);
                [MRE_b0, Var_b0] = BirthDeath_LNA(Theta_b0, tspan, 2, merr, y0, v0);
                [MRE_0b, Var_0b] = BirthDeath_LNA(Theta_0b, tspan, 2, merr, y0, v0);

                FD_MRE(i,j,d,:)     = (MRE_ff - MRE_f0 - MRE_0f + 2*MREa - MRE_b0 - MRE_0b + MRE_bb)/(2*Theta(i)*Theta(j)*delta^2);
                FD_Var(i,j,d,:,:)     = (Var_ff - Var_f0 - Var_0f + 2*Vara - Var_b0 - Var_0b + Var_bb)/(2*Theta(i)*Theta(j)*delta^2);
            end        
        end
        fprintf('\n')
    end
    toc

end

%% plot the second order F.D. MRE 
figure(hSens2_MRE)

k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        hold all
        plot(tspan, squeeze(FD_MRE(i,j,:)), 'o', 'MarkerSize', 1)       
        k=k+1;
    end
end

legend({'LNA++','F.D.'})
saveas(gcf,'~/Projects/LNA++/paper/Sens2MREComparison.pdf')

%% Plot Second order F.D. Variance

figure('Name','BirthDeath F.D. Second Order Sensitivities of the Autocovariance')
k=1;
d=3;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        imagesc(tspan,tspan,squeeze(FD_Var(i,j,1,:,:)))%, squeeze(clims(i,j,:))')
          if i==1
            title(labels{j})
        end        
%         if j==1
%             ylabel({labels{i},'Time (a.u.)'}, 'FontWeight', 'bold');
%         end
        set(gca,'FontSize',14)
        colorbar
%         xlabel('Time (a.u.)','FontWeight','bold')
        k=k+1;
    end
end

saveas(gcf,'~/Projects/LNA++/paper/Sens2Var_FD.pdf')

