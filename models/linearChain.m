%% performance testing using a system of simple linear conversions

% A -> B
% B -> C
% C -> D ... etc

% N species 
% reaction rate: concentration of upstream species x constant
% [ birth, death, conversions....]
clear all
F = @(phi,t,Theta) [ Theta(1), Theta(2)*phi(1), phi(1:end-1).*Theta(3:length(phi)+1), phi(2:end).*Theta(length(phi)+2:end)];
syms t

%%
computeTimeMat  = zeros(1,7);
runTimeMat0      = zeros(7,3);
runTimeMat1      = zeros(7,3);
runTimeMat2      = zeros(7,3);

%% estimate compute and run time for models of with 3-10 species
for k=1:7
    N=2+k;
    disp(N)
    phi     = sym('phi',[1,N]);
    Theta   = sym('Theta',[1,2*N]); 
    Theta   = sym(Theta, 'positive'); % parameters are positive
    SS = S_gen_linearChain(N);
    FF = matlabFunction(F(phi,t,Theta), 'vars', {phi, t, Theta});

    modelName = ['chain' int2str(N)];
    
    tic
    generateLNAComponents(modelName, SS, FF, phi, Theta, 'NONE')
    computeTimeMat(k) = toc;
    npar = length(Theta);
    compileLNA(modelName, SS, npar); 
    addpath(modelName);
    
    disp('Running simulations...')
    tspan       = 0:0.1:2.5;

    parfor i=1:10
        Theta   = rand(1,2*N)*100; %[1 0.5 2 3 4 5];
        Y0      = rand(1,N)*1000;
        
        % no sensitivities
        tic
        [M,S]     = eval(['chain' int2str(N) '_LNA(Theta, tspan, 1:N, 0, Y0, zeros(1,N*(N+1)/2))']);
        runTimeMat0(k,i) = toc;        
        fprintf('.')

        % first order sens.
        [M1,S1,dM,dS] = eval(['chain' int2str(N) '_LNA(Theta, tspan, 1:N, 0, Y0, zeros(1,N*(N+1)/2))']);
%             chain3_LNA(Theta, tspan, 1:N, 0, Y0);
        runTimeMat1(k,i) = toc;        
        fprintf('.')
        
        % second order sens.
        [M2,S2,dM,dS,d2M,d2S] = eval(['chain' int2str(N) '_LNA(Theta, tspan, 1:N, 0, Y0, zeros(1,N*(N+1)/2))']);
%             chain3_LNA(Theta, tspan, 1:N, 0, Y0);
        runTimeMat2(k,i) = toc;        
        fprintf('.')
    end
    fprintf('\n')
end

%% save compute and run times
save('linearChainTimes', 'computeTimeMat', 'runTimeMat0', 'runTimeMat1', 'runTimeMat2')

%% plot the compute times
close all
plot(3:9, computeTimeMat, '.-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 30, 'LineWidth', 2)
set(gca, 'YScale', 'log', 'FontSize', 14)
xlabel('Network size', 'FontSize', 20)
ylabel('Compute time (s)', 'FontSize', 20)

saveas(gcf, '~/Documents/Promotion/Papers/LNA++/bioinfo01/Figures/linearChain_computeTime.eps')

%% plot the run times
close all
figure, hold on
% errorbar(3:9, median(runTimeMat0'), min(runTimeMat0'), max(runTimeMat0'), '.-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 30, 'LineWidth', 2)
% errorbar(3:9, median(runTimeMat1'), min(runTimeMat1'), max(runTimeMat1'), 'o-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 10, 'LineWidth', 2)
% errorbar(3:9, median(runTimeMat2'), min(runTimeMat2'), max(runTimeMat2'), 's-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 15, 'LineWidth', 2)

plot(3:9, median(runTimeMat0'), 'k.-', 'MarkerSize', 30, 'LineWidth', 2)
plot(3:9, median(runTimeMat1'), 'ko-', 'MarkerSize', 10, 'LineWidth', 2)
plot(3:9, median(runTimeMat2'), 'ks-', 'MarkerSize', 15, 'LineWidth', 2)


set(gca, 'YScale', 'log', 'FontSize', 14)
xlabel('Network size', 'FontSize', 20)
ylabel('Median Run time (s)', 'FontSize', 20)

legend('northwest', {'LNA', '+ 1st Order Sens.', '+ 2nd Order Sens.'}, 'FontSize', 14)

saveas(gcf, '~/Documents/Promotion/Papers/LNA++/bioinfo01/Figures/linearChain_runTime.eps')

%%
saveas(gcf, '~/Documents/Promotion/Papers/LNA++/bioinfo01/Figures/linearChain_computeTime.png')

%%
modelSize0 = @(n) (n + n.*n+1/2 + n.^2);
modelSize1 = @(n) (n + n.*n+1/2 + n.^2).*(1 + 3*n);
modelSize2 = @(n) (n + n.*n+1/2 + n.^2).*(1 + 3*n + 9*n.^2);

%% plot normalized run times
close all
figure, hold on

plot(3:9, median(runTimeMat0') ./ modelSize0(3:9), '.-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 30, 'LineWidth', 2)
plot(3:9, median(runTimeMat1') ./ modelSize1(3:9), 'o-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 10, 'LineWidth', 2)
plot(3:9, median(runTimeMat2') ./ modelSize2(3:9), 's-', 'Color', [0.4 0.4 0.4], 'MarkerSize', 15, 'LineWidth', 2)

set(gca, 'YScale', 'log', 'FontSize', 14, 'YLim', [0.0005, 0.01])
xlabel('Network size', 'FontSize', 20)
ylabel('Normalized Median Run time (s)', 'FontSize', 20)

legend('northwest', {'LNA', '+ 1st Order Sens.', '+ 2nd Order Sens.'}, 'FontSize', 14)
saveas(gcf, '~/Documents/Promotion/Papers/LNA++/bioinfo01/Figures/linearChain_runTime_normalized.eps')
