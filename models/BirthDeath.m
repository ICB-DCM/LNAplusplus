% Simple birth/death model
% (C) Justin Feigelman
% justin.feigelman@helmholtz-muenchen.de
% Helmholtz Zentrum Muenchen

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

addpath('../matlab')
generateLNA('BirthDeath/BirthDeath.xml', 'BirthDeath', 'BOTH')

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

set(gcf, 'PaperPosition', [0.25, 2.5, 6,6])
labels = {'k_m','k_p','g_m','g_p'};
for i=1:4
    subplot(2,2,i)
    plot(tspan, Sens_MRE(i,:), 'LineWidth', 3)
    ylabel(labels{i})
    xlabel('Time')
    set(gca,'FontSize',20)    
end

%% First order sensitivity of the variance
hSens1_Var = figure('Name','BirthDeath First Order Sensitivities of the Autocovariance');
set(gcf, 'PaperPosition', [0.25, 2.5, 8,6])
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan,Sens_Var(:,:,i))
    title(['Sensitivity ' labels{i}])
%     xlabel('Time (a.u.)')
%     ylabel('Time (a.u.)')
    set(gca,'FontSize',20)
    h = colorbar;
    h.Label.String = ['Sensitivity ' labels{i}];
end
% colorbar
%%
saveas(gcf, '~/Projects/LNA++/paper/FirstOrderVarSensitivityAnalytical.pdf')


%% compute finite difference approximation to first order sensitivities
% delta=0.0001;
delta=5.9948e-04
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
    h = colorbar;
    h.Label.String = ['Sensitivity ' labels{i}];
end

saveas(gcf, '~/Projects/LNA++/paper/FirstOrderVarSensitivityFD.pdf')

%% subtract the variance sensitivities (normalized)
figure('Name', 'BirthDeath First Order Sensitivities Comparison')
set(gcf, 'PaperPosition', [0.25, 2.5, 8,6])
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan, abs( (squeeze(FD_Var(i,:,:)) - squeeze(Sens_Var(:,:,i))) ./ squeeze(Sens_Var(:,:,i))), [0,1])
    title(['Sensitivity ' labels{i}])
%     xlabel('Time (a.u.)')
%     ylabel('Time (a.u.)')
    set(gca,'FontSize',20)
    colorbar
end

%% subtract the variance sensitivities (raw)
figure('Name', 'BirthDeath First Order Sensitivities Comparison')
set(gcf, 'PaperPosition', [0.25, 2.5, 8,6])
for i=1:4
    subplot(2,2,i)
    imagesc(tspan,tspan, squeeze(FD_Var(i,:,:)) - squeeze(Sens_Var(:,:,i)))
    title(['Sensitivity ' labels{i}])
%     xlabel('Time (a.u.)')
%     ylabel('Time (a.u.)')
    set(gca,'FontSize',20)
    colorbar
end

%% try different perturbation sizes, compute the error
delta=logspace(-5,-1,10);
clear FD_MRE FD_Var var_err

for k=1:10
    for i=1:4
        [MREa, Vara] = BirthDeath_LNA(Theta, tspan, 2,0, y0, v0);
        Thetab = Theta;
        Thetab(i) = Thetab(i)*(1+delta(k));
        [MREb, Varb] = BirthDeath_LNA(Thetab, tspan, 2,0, y0, v0);    

%         FD_MRE(i,:)     = (MREb-MREa)/(Theta(i)*delta);
        FD_Var(i,:,:)   = (Varb-Vara)/(Theta(i)*delta(k));
        var_err(i,k) = norm( (squeeze(FD_Var(i,:,:)) - squeeze(Sens_Var(:,:,i))));
        fprintf('.')
    end
end

%% plot errors
figure
plot(var_err', 'o-','LineWidth',3)
set(gca, 'YScale', 'log')
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
        plot(tspan, squeeze(Sens2_MRE(i,j,:)), 'LineWidth', 3)

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
%%
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
 
% deltaVec = logspace(-2,-1,10);
% produces reasonable results with this perturbation size
% deltaVec=10^-1.3;
deltaVec=10^-4
for d=1:length(deltaVec)
    delta=deltaVec(d);

    tic
    for i=1:4
        for j=1:4
            fprintf('.')
            if i==j
                Theta_f = Theta; Theta_b=Theta;
%                 % relative perturbation
%                 Theta_f(i) = Theta(i)*(1+delta);
%                 Theta_b(i) = Theta(i)*(1-delta);
                
                % absolute peturbation
                Theta_f(i) = Theta(i)+delta;
                Theta_b(i) = Theta(i)-delta;

                [MRE_f, Var_f] = BirthDeath_LNA(Theta_f, tspan, 2, merr, y0, v0);
                [MRE_b, Var_b] = BirthDeath_LNA(Theta_b, tspan, 2, merr, y0, v0);    

%                 FD_MRE(i,i,d,:)     = (MRE_f+MRE_b-2*MREa)/(Theta(i)*delta)^2;
%                 FD_Var(i,i,d,:,:)   = (Var_f+Var_b-2*Vara)/(Theta(i)*delta)^2;
                
                FD_MRE(i,i,d,:)     = (MRE_f+MRE_b-2*MREa)/(delta^2);
                FD_Var(i,i,d,:,:)   = (Var_f+Var_b-2*Vara)/(delta^2);
            else
                Theta_ff = Theta; Theta_f0 = Theta; Theta_0f = Theta; 
                Theta_b0 = Theta; Theta_0b = Theta; Theta_bb = Theta;

                Theta_fb = Theta; Theta_bf = Theta;
                
%                 % relative perturbation size
%                 Theta_ff([i j]) = Theta([i j])*(1+delta);           
%                 Theta_f0(i) = Theta(i)*(1+delta);
%                 Theta_0f(j) = Theta(j)*(1+delta);
% 
%                 Theta_b0(i) = Theta(i)*(1-delta);
%                 Theta_0b(j) = Theta(j)*(1-delta);
%                 Theta_bb([i j]) = Theta([i j])*(1-delta);
% 
%                 Theta_fb([i j]) = Theta([i j]).*[1+delta,1-delta];
%                 Theta_bf([i j]) = Theta([i j]).*[1-delta,1+delta];
                
                % absolute perturbation size
                Theta_ff([i j]) = Theta([i j])+delta;           
                Theta_f0(i) = Theta(i)+delta;
                Theta_0f(j) = Theta(j)+delta;

                Theta_b0(i) = Theta(i)-delta;
                Theta_0b(j) = Theta(j)-delta;
                Theta_bb([i j]) = Theta([i j])-delta;

                Theta_fb([i j]) = Theta([i j])+ delta*[1,-1];
                Theta_bf([i j]) = Theta([i j])+ delta*[-1,1];
                
                [MRE_ff, Var_ff] = BirthDeath_LNA(Theta_ff, tspan, 2, merr, y0, v0);
%                 [MRE_f0, Var_f0] = BirthDeath_LNA(Theta_f0, tspan, 2, merr, y0, v0);
%                 [MRE_0f, Var_0f] = BirthDeath_LNA(Theta_0f, tspan, 2, merr, y0, v0);                
% 
                [MRE_bb, Var_bb] = BirthDeath_LNA(Theta_bb, tspan, 2, merr, y0, v0);
%                 [MRE_b0, Var_b0] = BirthDeath_LNA(Theta_b0, tspan, 2, merr, y0, v0);
%                 [MRE_0b, Var_0b] = BirthDeath_LNA(Theta_0b, tspan, 2, merr, y0, v0);

                [MRE_fb, Var_fb] = BirthDeath_LNA(Theta_fb, tspan, 2, merr, y0, v0);
                [MRE_bf, Var_bf] = BirthDeath_LNA(Theta_bf, tspan, 2, merr, y0, v0);
                
%                 FD_MRE(i,j,d,:)     = (MRE_ff - MRE_f0 - MRE_0f + 2*MREa - MRE_b0 - MRE_0b + MRE_bb)/(2*Theta(i)*Theta(j)*delta^2);
%                 FD_Var(i,j,d,:,:)   = (Var_ff - Var_f0 - Var_0f + 2*Vara - Var_b0 - Var_0b + Var_bb)/(2*Theta(i)*Theta(j)*delta^2);
                
%               alternative:
%                % relative perturbation
%                 FD_MRE(i,j,d,:,:)   = (MRE_ff - MRE_fb - MRE_bf + MRE_bb)/(4*Theta(i)*Theta(j)*delta^2);
%                 FD_Var(i,j,d,:,:)   = (Var_ff - Var_fb - Var_bf + Var_bb)/(4*Theta(i)*Theta(j)*delta^2);

                % using absolute perturbation
                FD_MRE(i,j,d,:,:)   = (MRE_ff - MRE_fb - MRE_bf + MRE_bb)/(4*delta^2);
                FD_Var(i,j,d,:,:)   = (Var_ff - Var_fb - Var_bf + Var_bb)/(4*delta^2);        
                
                
                var_err2(d,i,j) = norm( squeeze(FD_Var(i,j,d,:,:)) - squeeze(Sens2_Var(i,j,:,:)) ); 
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
        plot(tspan, squeeze(FD_MRE(i,j,:)), 'o', 'MarkerSize', 2)       
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

%%
saveas(gcf,'~/Projects/LNA++/paper/Sens2Var_FD.pdf')

%% comparison of variance sensitivities
figure('Name','BirthDeath Second Order Sensitivities comparison (normalized)')
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        imagesc(tspan,tspan, abs((squeeze(FD_Var(i,j,1,:,:)) - squeeze(Sens2_Var(i,j,:,:))) ./ squeeze(Sens2_Var(i,j,:,:))), [0,1] )
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

%% comparison of variance sensitivities unnormalized
figure('Name','BirthDeath Second Order Sensitivities comparison (raw)')
k=1;
for i=1:4
    for j=1:4
        subplot(4,4,k)
        imagesc(tspan,tspan, (squeeze(FD_Var(i,j,1,:,:)) - squeeze(Sens2_Var(i,j,:,:))))%  ./ squeeze(Sens2_Var(i,j,:,:)) )
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



%% show error of FD estimator for variance sens.

X = reshape(var_err2, 10, 16);

figure
plot(X, 'o-', 'LineWidth', 3)
set(gca, 'YScale', 'log')

%% compare FD error for different deltas
figure
for k=1:10
    subplot(3,4,k)    
    imagesc( squeeze(FD_Var(4,4,k,:,:)) - squeeze(Sens2_Var(4,4,:,:)))
end