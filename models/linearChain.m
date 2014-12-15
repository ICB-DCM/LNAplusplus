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
runTimeMat      = zeros(7,3);

%%
for k=3
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
    
    for i=1:3
        Theta   = rand(1,2*N)*100; %[1 0.5 2 3 4 5];
        Y0      = rand(1,N)*1000;
        
        tspan       = 0:0.1:10;
        fprintf('.')
        tic
        [M,S,dM,dS] = eval(['chain' int2str(N) '_LNA(Theta, tspan, 1:N, 0, Y0, zeros(1,N*(N+1)/2))']);
%             chain3_LNA(Theta, tspan, 1:N, 0, Y0);
        runTimeMat(k,i) = toc;        
    end
    fprintf('\n')
end
%%
close all
figure
plot(M')

figure
plot(squeeze(dM(:,1,:))')

%%
