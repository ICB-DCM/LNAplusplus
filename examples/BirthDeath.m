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

%%  define symbolic variables
syms k_m k_p g_m g_p real
syms m p real

%% define reaction fluxes
f = @(phi,t,Theta) ...
    [ Theta(1), Theta(3)*phi(1), Theta(2)*phi(1), Theta(4)*phi(2)];

phi     = [m,p]; % state vector
Theta   = [k_m, k_p, g_m, g_p]; % parameter vector
npar 	= length(Theta); % 4 parameters

%% generate the C code for the components of the LNA model
generateLNAComponents(model, S, f, phi, Theta)

%% compile code
compileLNA(model, S, npar); 

%%
%%%%%%%%%%%%%%% run some simulations and show the results %%%%%%%%%%%%%%%     
% parameters in same order as above
Theta = [20, 25, 10, 1];

% observation times for the simulation
tspan = 0:0.1:10;

% add mex file folder to matlab path
addpath('mex/')

% compute solution assuming steady state and no measurement error
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan);
toc

% plot MRE
close all
figure
subplot(121)
plot(tspan,MRE(2,:)) % just plot protein
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')

% plot variance
subplot(122)
imagesc(squeeze(Var(2,2,:,:)))
xlabel('Time')
ylabel('Time')
title('Protein covariance matrix')

%% specify an observed variable
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2); % outputs protein (species 2) only
toc

%% specify measurement error (same for each species)
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2, 50); % mRNA and protein have measurement variance 50
toc

figure
subplot(121)
plot(tspan,MRE)
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, Var)
xlabel('Time')
ylabel('Time')
title('Protein covariance matrix')

%% specify measurement error (each separately)
close all
tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, [1,2], [50, 100]); % measurement error 50 and 100 for mRNA and protein, respectively
toc

figure
subplot(121)
plot(tspan,MRE(2,:))
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, squeeze(Var(2,2,:,:)))
xlabel('Time')
ylabel('Time')
title('Protein covariance matrix')

%% specify initial values, variances

y0 = [2, 200]; % inital values
v0 = [10, 1, 500]; % variances
% v0 = [variance(mRNA), covariance(mRNA,protein), variance(protein)]

tic
[MRE, Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
toc

figure
subplot(121)
plot(tspan,MRE)
xlabel('Time')
ylabel('Protein')
title('Protein trajectory')
subplot(122)
imagesc(tspan,tspan, Var)
xlabel('Time')
ylabel('Time')
title('Protein covariance matrix')

%% 1st order sensitivities
tic
[MRE, Var, Sens_MRE, Sens_Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
% [MRE, Var, Sens_MRE, Sens_Var] = BirthDeath_LNA(tspan, Theta, merr, 2);
toc
%% first order sensitivity plots: analytical
close all
figure(1)
labels = {'k_m','k_p','g_m','g_p'};
for i=1:4
    subplot(2,2,i)
    plot(tspan, Sens_MRE(i,:), 'LineWidth', 2)
    title(['Sensitivity ' labels{i}])
    xlabel('Time')
    set(gca,'FontSize',20)    
end

figure(2)
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan,Sens_Var(:,:,i))
    title(['Sensitivity ' labels{i}])
    xlabel('Time')
    ylabel('Time')
    set(gca,'FontSize',20)
end

%%
%saveas(gcf, 'FirstOrderVarSensitivityAnalytical.pdf')


%% first order senstivity plots: finite difference approximation
delta=0.01;
clear FD_MRE FD_Var

for i=1:4
    [MREa, Vara] = BirthDeath_LNA(tspan, Theta, merr, 2, y0, v0);
    Thetab = Theta;
    Thetab(i) = Thetab(i)*(1+delta);
    [MREb, Varb] = BirthDeath_LNA(tspan, Thetab, merr, 2, y0, v0);    
    
    FD_MRE(i,:)     = (MREb-MREa)/(Theta(i)*delta);
    FD_Var(i,:,:)   = (Varb-Vara)/(Theta(i)*delta);
end


%% comparison of the MRE first order sensitivities
figure(1)
for i=1:4
    subplot(2,2,i)
    hold all
    plot(tspan, FD_MRE(i,:), 'o')
%     title(['Sensitivity ' labels{i}])
    xlabel('Time')
    set(gca,'FontSize',20)
end

legend('NorthEast', {'Analytical','Finite Difference'})
% saveas(gcf, 'FirstOrderMRESensitivityComparison.pdf')

%% plot finite difference first order var sensitivity

figure
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan,squeeze(FD_Var(i,:,:)))
    title(['Sensitivity ' labels{i}])
    xlabel('Time')
    ylabel('Time')
    set(gca,'FontSize',20)
end
%saveas(gcf, 'FirstOrderVarSensitivityFD.pdf')

%% second order sensitivities

tic
[MRE, Var, Sens_MRE, Sens_Var, Sens2_MRE, Sens2_Var] = BirthDeath_LNA(Theta, tspan, 2, 0, y0, v0);
toc
%% second order sensitivities assuming steady state
tic
[MRE, Var, Sens_MRE, Sens_Var, Sens2_MRE_SS, Sens2_Var_SS] = BirthDeath_LNA(Theta, tspan, 2, 0);
toc
%% plot second order sensitivity, MRE
close all
figure
labels = {'k_m','k_p','g_m','g_p','m_0','p_0'};
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        plot(tspan, squeeze(Sens2_MRE(i,j,:)), 'LineWidth', 2)

        xlabel('Time (a.u.)')
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
%% plot second order sensitivity, variance
close all
figure
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        imagesc(tspan,tspan,squeeze(Sens2_Var(i,j,:,:)))
%         title(sprintf('Sensitivity(%s,%s)', labels{i}, labels{j}))
        if i==1
            title(labels{j})
        end        
        if j==1
            ylabel({labels{i},'Time (a.u.)'}, 'FontWeight', 'bold');
        end
        set(gca,'FontSize',14)

%         ylabel('Time (a.u.)')
        xlabel('Time (a.u.)','FontWeight','bold')
        k=k+1;
    end
end

%% compute second order finite difference 
clear FD_MRE FD_Var
% Theta = [10 10 2 1];
% y0 = [5 50];
Theta = [20 25 10 1];
y0 = [2 200];

[MREa, Vara] = BirthDeath_LNA(tspan, Theta, merr, 2, y0, v0);

deltaVec = 10^-1.3; % produces reasonable results with this perturbation size
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

                [MRE_f, Var_f] = BirthDeath_LNA(tspan, Theta_f, merr, 2, y0, v0);
                [MRE_b, Var_b] = BirthDeath_LNA(tspan, Theta_b, merr, 2, y0, v0);    

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

                [MRE_ff, Var_ff] = BirthDeath_LNA(tspan, Theta_ff, merr, 2, y0, v0);
                [MRE_f0, Var_f0] = BirthDeath_LNA(tspan, Theta_f0, merr, 2, y0, v0);
                [MRE_0f, Var_0f] = BirthDeath_LNA(tspan, Theta_0f, merr, 2, y0, v0);                

                [MRE_bb, Var_bb] = BirthDeath_LNA(tspan, Theta_bb, merr, 2, y0, v0);
                [MRE_b0, Var_b0] = BirthDeath_LNA(tspan, Theta_b0, merr, 2, y0, v0);
                [MRE_0b, Var_0b] = BirthDeath_LNA(tspan, Theta_0b, merr, 2, y0, v0);

                FD_MRE(i,j,d,:)     = (MRE_ff - MRE_f0 - MRE_0f + 2*MREa - MRE_b0 - MRE_0b + MRE_bb)/(2*Theta(i)*Theta(j)*delta^2);
                FD_Var(i,j,d,:,:)     = (Var_ff - Var_f0 - Var_0f + 2*Vara - Var_b0 - Var_0b + Var_bb)/(2*Theta(i)*Theta(j)*delta^2);
            end        
        end
        fprintf('\n')
    end
    toc

end

%% plot the second order F.D. MRE 
figure(1)

k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        hold all
        plot(tspan, squeeze(FD_MRE(i,j,:)), 'o', 'MarkerSize', 3)       
        k=k+1;
    end
end

legend({'Analytical','Finite Difference'})
%saveas(gcf,'Sens2MREComparison.pdf')

%% Plot Second order F.D. Variance

figure
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        imagesc(tspan,tspan,squeeze(FD_Var(i,j,1,:,:)))
          if i==1
            title(labels{j})
        end        
        if j==1
            ylabel({labels{i},'Time (a.u.)'}, 'FontWeight', 'bold');
        end
        set(gca,'FontSize',14)

        xlabel('Time (a.u.)','FontWeight','bold')
        k=k+1;
    end
end

%%



