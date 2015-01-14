function output = generateLNAComponents(modelName, S, F, phi, Theta, varargin) 
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
% optional arguments:
%  computeSS: Can be 'Both' (default), 'Y0', 'V0', or 'None'.  
%   'Y0': attempt to compute the steady state solution to the MRE
%   'V0': attempt to compute the steady state solution to the Variance
%   'Both': compute both (default)
%   'None': do not attempt to compute the steady states.

%% parse input
if nargin > 5
    computeSS = varargin{1};
    switch computeSS
        case 'Y0'
            COMPUTE_Y0 = true;
            COMPUTE_V0 = false;
        case 'V0'
            COMPUTE_Y0 = false;
            COMPUTE_V0 = true;
        case 'BOTH'
            COMPUTE_Y0 = true;
            COMPUTE_V0 = true;
        case 'NONE'
            COMPUTE_Y0 = false;
            COMPUTE_V0 = false;
        otherwise
            error('computeSS must be either ''Y0'', ''V0'', ''BOTH'' or ''NONE''')
    end
else
    COMPUTE_V0 = false;
    COMPUTE_Y0 = false;        
end
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

reactionFlux = F(phi, t, Theta);

% jacobian of flux vector
J           = Jacobian(reactionFlux, phi);

% df/dtheta
dFdTheta    = Jacobian(reactionFlux, Theta);

% d2f/dtheta2
d2fdTheta2  = Jacobian(dFdTheta, Theta);

% A
Afunc           = S*J;

% sensitivities of A
dAdTheta    = Jacobian(Afunc, Theta);
dAdPhi      = Jacobian(Afunc,phi);
d2AdTheta2  = Jacobian(dAdTheta, Theta);
d2AdPhi2    = Jacobian(dAdPhi, phi);

% E
Efunc           = S*sqrt(diag(reactionFlux));

% sensitivities of E
dEdTheta    = Jacobian(Efunc, Theta);
d2EdTheta2   = Jacobian(dEdTheta, Theta);
dEdPhi      = Jacobian(Efunc, phi);
d2EdPhi2    = Jacobian(dEdPhi, phi);

%% solve for the initial steady state 
nvar = length(phi);
npar = length(Theta)+nvar;

if COMPUTE_Y0
    Y0 = solveSS_mre(S,reactionFlux,phi);
    if isempty(Y0)
        error('Could not compute symbolic steady state.  Please supply initial condition');
    end
else
    Y0 = zeros(1,nvar); % placeholder
end
% [V0 systemJacobian MI] = solveSS_var(A,E,F,S,phi,Theta,sym(Y0));
% [V0 systemJacobian] = solveSS_var(A,E,F,S,phi,Theta,Y0);

if COMPUTE_V0
    V0 = solveSS_var(Afunc,Efunc,reactionFlux,S,phi,Theta,phi0);
    if isempty(V0)
        error('Could not compute symbolic steady state.  Please supply initial condition');
    end
else
    V0 = zeros(1, nvar*(nvar+1)/2); % placeholder
end

%% system Jacobian
sysVar = phi;
V = sym('V', [nvar nvar]);
V = sym(V, 'real');
for i=1:length(phi)
    for j=1:i
        sysVar = [sysVar, V(i,j)];
    end
end
% fundamental matrix
Phi = sym('Phi', [numel(phi), numel(phi)]);
Phi = sym(Phi,'real');
sysVar = [sysVar reshape(Phi,1,[])];
dVdt = Afunc*V + V*Afunc' + subs(Efunc*Efunc.');
RHS = [S*reactionFlux'; dVdt(find(triu(ones(size(dVdt))))); reshape(Afunc*Phi,[],1)];

systemJacobian = Jacobian(RHS, sysVar);
systemJacobian_diag = diag(systemJacobian); % just the diagonal

% preconditioner
% gamma=sym('gamma','real');
% M=eye(size(systemJacobian))-gamma*systemJacobian;

% MI = pinv(M);

%% initial sensitivities
disp('computing initial sensitivities')

if COMPUTE_Y0
    % dY/dTheta
    S0  = Jacobian(Y0, Theta);
    % d2Y/dTheta2
    S20 = Jacobian(S0, Theta);
else
   S0   = zeros(1,nvar*npar); % placeholder
   S20  = zeros(1,nvar*npar*npar); % placeholder
end

if COMPUTE_V0
    % dV/dTheta
    SV0 = Jacobian(V0, Theta);
    % d2V/dTheta2
    S2V0 = Jacobian(SV0, Theta);
else
    SV0     = zeros(1,nvar*(nvar+1)/2*npar); % placeholder
    S2V0    = zeros(1,nvar*(nvar+1)/2*npar*npar); % placeholder
end
%% mixed second derivatives
disp('Computing Mixed Second Derivatives')
d2AdPhidTheta = Jacobian(dAdPhi, Theta);
d2AdThetadPhi = Jacobian(dAdTheta, phi);

d2EdPhidTheta = Jacobian(dEdPhi, Theta);
d2EdThetadPhi = Jacobian(dEdTheta, phi);

%% generate the matlab functions based on these symbolic functions
% disp('Generating m-files')
% clear t, t=sym('t','real');
% matlabFunction(reactionFlux, 'file', [dirName '/matlab/reactionFlux'], 'vars', {phi, t, Theta});
% matlabFunction(J, 'file', [dirName '/matlab/J'], 'vars', {phi, t, Theta});
% matlabFunction(dFdTheta, 'file', [dirName '/matlab/dFdTheta'], 'vars', {phi, t, Theta});
% matlabFunction(d2fdTheta2, 'file', [dirName '/matlab/d2fdTheta2'], 'vars', {phi, t, Theta});
% matlabFunction(Afunc, 'file', [dirName '/matlab/Afunc'], 'vars', {phi, t, Theta});
% matlabFunction(dAdTheta, 'file', [dirName '/matlab/dAdTheta'], 'vars', {phi, t, Theta});
% matlabFunction(dAdPhi, 'file', [dirName '/matlab/dAdPhi'], 'vars', {phi, t, Theta});
% matlabFunction(d2AdTheta2, 'file', [dirName '/matlab/d2AdTheta2'], 'vars', {phi, t, Theta});
% matlabFunction(d2AdPhi2, 'file', [dirName '/matlab/d2AdPhi2'], 'vars', {phi, t, Theta});
% matlabFunction(Efunc, 'file', [dirName '/matlab/Efunc'], 'vars', {phi, t, Theta});
% matlabFunction(dEdTheta, 'file', [dirName '/matlab/dEdTheta'], 'vars', {phi, t, Theta});
% matlabFunction(d2EdTheta2, 'file', [dirName '/matlab/d2EdTheta2'], 'vars', {phi, t, Theta});
% matlabFunction(dEdPhi, 'file', [dirName '/matlab/dEdPhi'], 'vars', {phi, t, Theta});
% matlabFunction(d2EdPhi2, 'file', [dirName '/matlab/d2EdPhi2'], 'vars', {phi, t, Theta});
% matlabFunction(systemJacobian, 'file', [dirName '/matlab/systemJacobian'], 'vars', {phi, t, Theta});
% matlabFunction(systemJacobian_diag, 'file', [dirName '/matlab/systemJacobian_diag'], 'vars', {phi, t, Theta});
% 
% % gamma=sym('gamma','real');
% % matlabFunction(MI, 'file', [dirName '/matlab/MI'], 'vars', {phi, t, Theta, gamma});
% 
% matlabFunction(sym(Y0), 'file', [dirName '/matlab/Y0'], 'vars', {Theta});
% matlabFunction(sym(V0), 'file', [dirName '/matlab/V0'], 'vars', {Theta});
% 
% matlabFunction(sym(S0), 'file', [dirName '/matlab/S0'], 'vars', {Theta});
% matlabFunction(sym(SV0), 'file', [dirName '/matlab/SV0'], 'vars', {Theta});
% matlabFunction(sym(S20), 'file', [dirName '/matlab/S20'], 'vars', {Theta});
% matlabFunction(sym(S2V0), 'file', [dirName '/matlab/S2V0'], 'vars', {Theta});
% 
% matlabFunction(d2EdPhidTheta, 'file', [dirName '/matlab/d2EdPhidTheta'], 'vars', {phi, t, Theta});
% matlabFunction(d2EdThetadPhi, 'file', [dirName '/matlab/d2EdThetadPhi'], 'vars', {phi, t, Theta});
% 
% matlabFunction(d2AdPhidTheta, 'file', [dirName '/matlab/d2AdPhidTheta'], 'vars', {phi, t, Theta});
% matlabFunction(d2AdThetadPhi, 'file', [dirName '/matlab/d2AdThetadPhi'], 'vars', {phi, t, Theta});

%% generate the C library
disp('Generating C source')
% add the codegen tags to all the generated m files
% codegenify( [dirName '/matlab'] );

% % TODO: use the codegen configuration object to dynamically specify
% % dimensions of input arguments phi and theta
olddir = cd([dirName '/matlab']);

%     ' d2EdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...

% codegen_cmd = sprintf(['codegen -v -O disable:openmp -c -config:lib -d %s '...
%     ' reactionFlux -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' J -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' dFdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2fdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' Afunc -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' dAdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' dAdPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2AdPhi2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    
%     ' d2AdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2AdThetadPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2AdPhidTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' Efunc -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' dEdTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2EdTheta2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2EdThetadPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    
%     ' d2EdPhidTheta -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ... 
%     ' dEdPhi -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' d2EdPhi2 -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...    ' MI -args {zeros(1,NVAR),0,zeros(1,NPAR),0}' ...     
%     ' S0 -args {zeros(1,NPAR)}' ...
%     ' S20 -args {zeros(1,NPAR)}' ...
%     ' SV0 -args {zeros(1,NPAR)}' ...
%     ' S2V0 -args {zeros(1,NPAR)}' ...      
%     ' Y0 -args {zeros(1,NPAR)}' ...    
%     ' V0 -args {zeros(1,NPAR)}' ... 
%     ' systemJacobian_diag -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' systemJacobian -args {zeros(1,NVAR),0,zeros(1,NPAR)}' ...
%     ' -I /usr/include/c++/4.2.1/ -I /usr/include'], ...
%     '../C');
% 
% codegen_cmd = strrep(codegen_cmd, 'NVAR', int2str(length(phi)));
% codegen_cmd = strrep(codegen_cmd, 'NPAR', int2str(length(Theta)));
% 
% eval(codegen_cmd)
% 
% NVAR = int2str(length(phi));
% NPAR = int2str(length(Theta));
% 
% objs1 = {'reactionFlux', 'J', 'dFdTheta', 'd2fdTheta2', 'Afunc', 'dAdTheta', ...
%     'dAdPhi', 'd2AdPhi2', 'd2AdTheta2', 'd2AdThetadPhi', 'd2AdPhidTheta', ...
%     'd2AdPhidTheta', 'Efunc', 'dEdTheta', 'd2EdTheta2', 'd2EdThetadPhi', ...
%     'd2EdPhidTheta', 'dEdPhi', 'd2EdPhi2', 'systemJacobian', 'systemJacobian_diag'};
% 
% for o=objs1
%     fprintf('%s ', o{:})
%     tic
%     eval(sprintf('codegen -c -d ../C %s -args {zeros(1,%s),0,zeros(1,%s)}', o{:}, NVAR, NPAR))
%     t=toc;
%     fprintf('%0.3f\n', t);
% end
% 
% objs2 = {'S0','S20','SV0','S2V0','Y0','V0'};
% for o = objs2
%     fprintf('%s ', o{:});
%     tic
%     eval(sprintf('codegen -c -d ../C %s -args {zeros(1,%s)}', o{:}, NPAR))
%     t=toc;
%     fprintf('%0.3f\n', t);
% end

%% generate C code

% requires phi, t, Theta
objs1 = {'reactionFlux', 'J', 'dFdTheta', 'd2fdTheta2', 'Afunc', 'dAdTheta', ...
    'dAdPhi', 'd2AdPhi2', 'd2AdTheta2', 'd2AdThetadPhi', 'd2AdPhidTheta', ...
    'd2AdPhidTheta', 'Efunc', 'dEdTheta', 'd2EdTheta2', 'd2EdThetadPhi', ...
    'd2EdPhidTheta', 'dEdPhi', 'd2EdPhi2', 'systemJacobian', 'systemJacobian_diag'};

% requires only Theta
objs2 = {'S0','S20','SV0','S2V0','Y0','V0'};
addpath('../../../matlab')
for i = 1:length(objs1)
    fprintf('Generating %s source...\n', objs1{i})
    genCCode(eval(objs1{i}), objs1{i}, {phi, t, Theta});
end

for i = 1:length(objs2)
    fprintf('Generating %s source...\n', objs2{i})
    genCCode(sym(eval(objs2{i})), objs2{i}, {Theta});  % make sure the zero matrices are symbolic vars
end    


%% generate MODEL_DEF file
disp('Generating model definition file')
f = fopen('../C/MODEL_DEF.h', 'w');
fprintf(f, '#define STOICH_MAT %s\n', strjoin(cellfun(@num2str,num2cell(reshape(S',1,[])), 'UniformOutput', false),','));
fprintf(f, '#define NVAR %d\n', length(phi));
fprintf(f, '#define NPAR %d\n', length(Theta)-length(phi));
fprintf(f, '#define NREACT %d\n', size(S,2));
if COMPUTE_Y0
    fprintf(f, '#define COMPUTE_Y0\n');
end
if COMPUTE_V0
    fprintf(f, '#define COMPUTE_V0\n');
end

fclose(f);

cd(olddir)

disp('Done')

end

function Y0=solveSS_mre(S,F,phi) %#codegen
% solve for the steady state with the parameters specified

F = reshape(F,[],1);
% tmp = solve(S*F,phi);
tmp1=num2cell(S*F);
tmp2=num2cell(phi);
tmp = solve(tmp1{:},tmp2{:}); % solve for phi

if isempty(tmp)
    warning('Could not solve for steady state initial conditions.  Please supply the initial condition explicitly.')
    Y0 = [];
    return
end
tmp2 = struct2cell(tmp);
Y0 = [tmp2{:}];

% for i=1:length(phi0), Y0(i)=tmp.(char(phi0(i))); end

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

f1 = matlabFunction(dVdt, 'vars', {phi, t, Theta, V}); % function handle to dVdt
f2 = f1(Y0,t,Theta,V); % matrix of symbols
f3 = f2(find(triu(ones(size(f2)))));
% convert to cell arrays for backwards compatibility
tmp1=num2cell(f3);
tmp2=num2cell(V(find(triu(V))));
V0 = solve(tmp1{:},tmp2{:});

if isempty(V0)
    warning('Could not solve for steady state variance.  Setting initial variance to zero.')
    V0 = zeros(nvar*(nvar+1)/2,1);
else
    V0 = subs(V(find(triu(V))), V0); 
end

% convert to symbolic vector


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
