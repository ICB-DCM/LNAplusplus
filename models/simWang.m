%%

function simWang()

%% REACTIONS

% R1: DNA + P2 ? DNA.P2
% R2: DNA.P2 ? DNA + P2
% R3: DNA ? DNA + P 
% R4: RNA ?> 0
% R5: 2 P ? P2
% R6: P2 ? 2 P
% R7: RNA ? RNA + P 
% R8: P ? 0
addpath('../matlab')
model = 'Wang';

syms DNA DNAPP M P PP
k = sym('k',[1,8]);

phi = [DNA, DNAPP, M, P, PP];
sbmlModel = TranslateSBML('Wang/Wang2010.xml');
[ S, ~, P ] = SBML2StoichProp(sbmlModel, true)

f = matlabFunction(P, 'Vars', {'phi','t','Theta'});

% S = [ -1, 1, 0, 0, -1;
%     1, -1, 0, 0, 1;
%     0, 0, 1, 0, 0;
%     0, 0, -1, 0, 0;
%     0, 0, 0, -2, 1;
%     0, 0, 0, 2, -1;
%     0, 0, 0, 1, 0;
%     0, 0, 0, -1, 0]';

    
%% flux vector

% generate the C code for the components of the LNA model
generateLNAComponents(model, S, f, phi, k, 'NONE')
% compile code
compileLNA(model); 

% function f = flux(phi, t, Theta)
% DNA = phi(1);
% DNAPP = phi(2);
% M = phi(3);
% P = phi(4);
% PP = phi(5);
% 
% f = [DNA*PP, DNAPP, DNA, M, P^2, PP, M, P] .* Theta(1:8);
% 

