%% Wang model
% Theta = 10*rand(1,8);
Theta = [0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1];
tspan = linspace(0,150,10);
Y0 = [200, 0, 0, 0, 0];
V0 = toLinear(zeros(5));
%%
tic
[MRE,Var] = Wang_LNA(Theta, tspan, 1:5, 0, Y0, V0);
toc
%%
labels = {'DNA', 'DNAP2', 'RNA', 'P', 'P2'};
close all
for k=1:5
    figure
    title(labels{k})
    shadedErrorBar(tspan, MRE(k,:), sqrt(diag(squeeze(Var(k,k,:,:)))))
end

%%
for k=1:5
    figure(k)
    title(labels{k})
end