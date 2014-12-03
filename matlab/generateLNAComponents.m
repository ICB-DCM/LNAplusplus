function output = generateLNAComponents(modelName, S, reactionFlux, phi, Theta) 
% generate all of the necessary components for the LNA model for the model
% with the given stoichiometric matrix S, and reactino flux vector
% reactionFlux
% input:
%  modelName: name of the model (string)
%  S: stoichiometric matrix
%  reactionFlux: function handle for computing the reaction fluxes for the
%   given state and parameters
%  phi: symbolic vector containing the macroscopic state variables
%  Theta: the model parameters (symbolic)

%% create directories to store the outputted scripts
dirName = modelName;
if ~exist(dirName, 'dir')
    mkdir(dirName)
end
if ~exist([dirName '/matlab'], 'dir')
    mkdir([dirName '/matlab']) % generated matlab functions
end
if ~exist([dirName '/C'], 'dir')
    mkdir([dirName '/C']) % generated C code
    %mkdir([dirName '/mex'])
end

disp('Computing symbolic derivatives')
t = sym('t', 'real');

%% add initial conditions to the vector of model parameters
phi0 = sym('phi0', [1, length(phi)]);
phi0 = sym(phi0, 'real');

Theta = [Theta, phi0]; % necessary for sensitivities wrt initial conditions

F = reactionFlux(phi, t, Theta);

% jacobian of flux vector
J           = Jacobian(F, phi);

% df/dtheta
dFdTheta    = Jacobian(F, Theta);

% d2f/dtheta2
d2fdTheta2  = Jacobian(dFdTheta, Theta);

% A
A           = S*J;

% sensitivities of A
dAdTheta    = Jacobian(A, Theta);
dAdPhi      = Jacobian(A,phi);
d2AdTheta2  = Jacobian(dAdTheta, Theta);
d2AdPhi2    = Jacobian(dAdPhi, phi);

% E
E           = S*sqrt(diag(F));

% sensitivities of E
dEdTheta    = Jacobian(E, Theta);
d2EdTheta2   = Jacobian(dEdTheta, Theta);
dEdPhi      = Jacobian(E, phi);
d2EdPhi2    = Jacobian(dEdPhi, phi);

%% solve for the initial steady state 
Y0 = solveSS_mre(S,F,phi,phi0);
if isempty(Y0)
    Y0 = @(Theta) error('Could not compute symbolic steady state.  Please supply initial condition');
end
% [V0 systemJacobian MI] = solveSS_var(A,E,F,S,phi,Theta,sym(Y0));
% [V0 systemJacobian] = solveSS_var(A,E,F,S,phi,Theta,Y0);
V0 = solveSS_var(A,E,F,S,phi,Theta,Y0);
%% initial sensitivities
disp('computing initial sensitivities')

% dY/dTheta
S0  = Jacobian(Y0, Theta);
% d2Y/dTheta2
S20 = Jacobian(S0, Theta);

% dV/dTheta
SV0 = Jacobian(V0, Theta);
% d2V/dTheta2
S2V0 = Jacobian(SV0, Theta);

%% mixed second derivatives

d2AdPhidTheta = Jacobian(dAdPhi, Theta);
d2AdThetadPhi = Jacobian(dAdTheta, phi);

d2EdPhidTheta = Jacobian(dEdPhi, Theta);
d2EdThetadPhi = Jacobian(dEdTheta, phi);

%% generate the matlab functions based on these symbolic functions
disp('Generating m-files')
clear t, t=sym('t','real');
matlabFunction(F, 'file', [dirName '/matlab/reactionFlux'], 'vars', {phi, t, Theta});
matlabFunction(J, 'file', [dirName '/matlab/J'], 'vars', {phi, t, Theta});
matlabFunction(dFdTheta, 'file', [dirName '/matlab/dFdTheta'], 'vars', {phi, t, Theta});
matlabFunction(d2fdTheta2, 'file', [dirName '/matlab/d2fdTheta2'], 'vars', {phi, t, Theta});
matlabFunction(A, 'file', [dirName '/matlab/Afunc'], 'vars', {phi, t, Theta});
matlabFunction(dAdTheta, 'file', [dirName '/matlab/dAdTheta'], 'vars', {phi, t, Theta});
matlabFunction(dAdPhi, 'file', [dirName '/matlab/dAdPhi'], 'vars', {phi, t, Theta});
matlabFunction(d2AdTheta2, 'file', [dirName '/matlab/d2AdTheta2'], 'vars', {phi, t, Theta});
matlabFunction(d2AdPhi2, 'file', [dirName '/matlab/d2AdPhi2'], 'vars', {phi, t, Theta});
matlabFunction(E, 'file', [dirName '/matlab/Efunc'], 'vars', {phi, t, Theta});
matlabFunction(dEdTheta, 'file', [dirName '/matlab/dEdTheta'], 'vars', {phi, t, Theta});
matlabFunction(d2EdTheta2, 'file', [dirName '/matlab/d2EdTheta2'], 'vars', {phi, t, Theta});
matlabFunction(dEdPhi, 'file', [dirName '/matlab/dEdPhi'], 'vars', {phi, t, Theta});
matlabFunction(d2EdPhi2, 'file', [dirName '/matlab/d2EdPhi2'], 'vars', {phi, t, Theta});

% matlabFunction(systemJacobian, 'file', [dirName '/matlab/systemJacobian'], 'vars', {phi, t, Theta});

gamma=sym('gamma','real');
% matlabFunction(MI, 'file', [dirName '/matlab/MI'], 'vars', {phi, t, Theta, gamma});

matlabFunction(sym(Y0), 'file', [dirName '/matlab/Y0'], 'vars', {Theta});
matlabFunction(sym(V0), 'file', [dirName '/matlab/V0'], 'vars', {Theta});

matlabFunction(sym(S0), 'file', [dirName '/matlab/S0'], 'vars', {Theta});
matlabFunction(sym(SV0), 'file', [dirName '/matlab/SV0'], 'vars', {Theta});
matlabFunction(sym(S20), 'file', [dirName '/matlab/S20'], 'vars', {Theta});
matlabFunction(sym(S2V0), 'file', [dirName '/matlab/S2V0'], 'vars', {Theta});

matlabFunction(d2EdPhidTheta, 'file', [dirName '/matlab/d2EdPhidTheta'], 'vars', {phi, t, Theta});
matlabFunction(d2EdThetadPhi, 'file', [dirName '/matlab/d2EdThetadPhi'], 'vars', {phi, t, Theta});

matlabFunction(d2AdPhidTheta, 'file', [dirName '/matlab/d2AdPhidTheta'], 'vars', {phi, t, Theta});
matlabFunction(d2AdThetadPhi, 'file', [dirName '/matlab/d2AdThetadPhi'], 'vars', {phi, t, Theta});

%% generate the C library
disp('Generating C source')
% add the codegen tags to all the generated m files
codegenify( [dirName '/matlab'] );

% TODO: use the codegen configuration object to dynamically specify
% dimensions of input arguments phi and theta
olddir = cd([dirName '/matlab']);

%     ' d2EdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...

codegen_cmd = sprintf(['codegen -config:lib -d %s '...
    ' reactionFlux -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' J -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' dFdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2fdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' Afunc -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' dAdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' dAdPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2AdPhi2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    
    ' d2AdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2AdThetadPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2AdPhidTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' Efunc -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' dEdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2EdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2EdThetadPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    
    ' d2EdPhidTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ... 
    ' dEdPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' d2EdPhi2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    ' MI -args {zeros(1,NVAR),0,zeros(1,NPAR),0}' ...     
    ' S0 -args {zeros(1,NPAR)}' ...
    ' S20 -args {zeros(1,NPAR)}' ...
    ' SV0 -args {zeros(1,NPAR)}' ...
    ' S2V0 -args {zeros(1,NPAR)}' ...      
    ' Y0 -args {zeros(1,NPAR)}' ...    
    ' V0 -args {zeros(1,NPAR)}' ... ' systemJacobian -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
    ' -I /usr/include/c++/4.2.1/ -I /usr/include'], ...
    '../C');
codegen_cmd = strrep(codegen_cmd, 'NVAR', int2str(length(phi)));
codegen_cmd = strrep(codegen_cmd, 'NPAR', int2str(length(Theta)));

eval(codegen_cmd);

%% generate MODEL_DEF file
f = fopen('../C/MODEL_DEF.h', 'w');
fprintf(f, '#define STOICH_MAT %s\n', strjoin(cellfun(@num2str,num2cell(reshape(S',1,[])), 'UniformOutput', false),','));
fprintf(f, '#define NVAR %d\n', length(phi));
fprintf(f, '#define NPAR %d\n', length(Theta)-length(phi));
fprintf(f, '#define NREACT %d\n', size(S,2));
fclose(f);

cd(olddir)

disp('Done')

end

function Y0=solveSS_mre(S,F,phi,phi0) %#codegen
% solve for the steady state with the parameters specified

t=0;
F = reshape(F,[],1);
% set the macroscopic flux to zero

fluxCellArray = mat2cell((S*subs(F, num2cell(phi), num2cell(phi0)))',1,ones(size(phi)));
phiCellArray  = mat2cell(phi0,1,ones(size(phi0)));

tmp=solve(fluxCellArray{:}, phiCellArray{:});

if isempty(tmp)
    warning('Could not solve for steady state initial conditions.  Please supply the initial condition explicitly.')
    Y0 = [];
    return
end

for i=1:length(phi0), Y0(i)=tmp.(char(phi0(i))); end

end

function [V0] = solveSS_var(A,E,F,S,phi,Theta,Y0) %#codegen
% solve for the steady state covariance matrix
% set the fluctuations to steady-state

% dVdt = AV + VA' + EE'
V = sym('V', [numel(phi), numel(phi)]);
V = sym(V, 'real');
V2 = V.';
V = subs(V2,V2(find(triu(V2,1))),V(find(triu(V,1)))); % symmetrize

dVdt = A*V + V*A' + subs(E*E.');

t = sym('t','real');
% set the upper triangular components to zero flux!

nvar = numel(phi);
%     VfluxCellArray  = mat2cell( reshape( dVdt, 1, []), 1, ones(1,nvar^2));
%     VarCellArray    = mat2cell( reshape( V, 1, []), 1, ones(1,nvar^2));
% 
%     V0              = solve(VfluxCellArray{:}, VarCellArray{:});

f1 = matlabFunction(dVdt, 'var', {phi, t, Theta, V}); % function handle to dVdt
f2 = f1(Y0,t,Theta,V); % matrix of symbols
f3 = f2(find(triu(ones(size(f2)))));
V0 = solve(f3, V(find(triu(V))));

if isempty(V0)
    warning('Could not solve for steady state variance.  Setting initial variance to zero.')
    V0 = zeros(nvar);
end

% convert to symbolic vector
V0 = subs(V(find(triu(V))), V0); 


% tmp = struct2cell(V0);
% V0              = [tmp{:}];
% V0              = reshape(V0,nvar,nvar);
% V0              = V0(find(triu(ones(size(V0))))); % just upper triangular region

% substitute the solution to the initial steady state of the species

%     varNames = cellfun(@(x) [char(x)], mat2cell(phi,1,ones(1,length(phi))), ...
%         'UniformOutput', false);
%     syms(varNames{:}, 'real'); % generate symbols in this scope
% 
%     thetaNames = cellfun(@(x) [char(x)], mat2cell(Theta,1,ones(1,length(Theta))), ...
%         'UniformOutput', false);
%     syms(thetaNames{:}, 'real'); % generate symbols in this scope
% 
%     for i=1:length(phi), sym(phi(i), 'real'), eval([char(phi(i)), '=', char(Y0(i))]), end
%     V0              = subs(V0);

% system Jacobian
% sysVar = phi;
% for i=1:length(phi)
%     for j=1:i
%         sysVar = [sysVar, V(i,j)];
%     end
% end
% Phi = sym('Phi', [numel(phi), numel(phi)]);
% Phi = sym(Phi,'real');
% 
% sysVar = [sysVar reshape(Phi,1,[])];
% 
% RHS = [S*F'; dVdt(find(triu(ones(size(dVdt))))); reshape(A*Phi,[],1)];
% 
% systemJacobian = Jacobian(RHS, sysVar);

% % preconditioner
% gamma=sym('gamma','real');
% M=eye(size(systemJacobian))-gamma*systemJacobian;
% 
% MI = pinv(M);

end

function J=Jacobian(x,y)

% J = jacobian(reshape(x,1,[]),y);
J = jacobian(reshape(x,[],1),y);
end
