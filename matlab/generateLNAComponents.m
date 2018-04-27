function output = generateLNAComponents(modelName, S, reactionFluxFun, phi, Theta, varargin) 
% generate all of the necessary components for the LNA model for the model
% with the given stoichiometric matrix S, and reactino flux vector
% reactionFlux
% input:
%  modelName: name of the model (string)
%  S: stoichiometric matrix
%  reactionFluxFun: function handle @(phi, t, Theta) for symbolically computing the reaction fluxes for the
%   given state and parameters
%  phi: symbolic vector containing the macroscopic state variables
%  Theta: symbolic vector containing the model parameters
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
lnaMatlabDir = fullfile(fileparts(mfilename('fullpath')));
dirName = fullfile(lnaMatlabDir, '..', 'models', modelName);
if ~exist(fullfile(fileparts(dirName)), 'dir')
    mkdir(fullfile(fileparts(dirName)))
end
if ~exist(dirName, 'dir')
    mkdir(dirName)
end
if ~exist([dirName '/matlab'], 'dir')
    mkdir([dirName '/matlab']) % generated matlab functions
end
if ~exist([dirName '/C'], 'dir')
    mkdir([dirName '/C']) % generated C code
end

disp('Computing symbolic derivatives')
t = sym('t', 'real');

%% add initial conditions to the vector of model parameters
phi0 = sym('phi0', [1, length(phi)]);
Theta = [Theta, phi0]; % necessary for sensitivities wrt initial conditions

assume(Theta>0)
assume(phi>0)
assume(phi0>0)

reactionFlux = reactionFluxFun(phi, t, Theta);

% jacobian of flux vector
J           = Jacobian(reactionFlux, phi);

% df/dtheta
dFdTheta    = Jacobian(reactionFlux, Theta);

% d2f/dtheta2
d2fdTheta2  = Jacobian(dFdTheta, Theta);

% A and sensitivities
Afunc          = S*J;
dAdTheta       = Jacobian(Afunc, Theta);
dAdPhi         = Jacobian(Afunc, phi);
d2AdTheta2     = Jacobian(dAdTheta, Theta);
d2AdPhi2       = Jacobian(dAdPhi, phi);
d2AdPhidTheta  = simplify(Jacobian(dAdPhi, Theta));
d2AdThetadPhi  = simplify(Jacobian(dAdTheta, phi));

% E*E^T and sensitivities
Efunc          = S*sqrt(diag(reactionFlux));
EEfunc         = simplify(Efunc*Efunc.');
dEEdTheta      = simplify(Jacobian(EEfunc, Theta));
dEEdPhi        = simplify(Jacobian(EEfunc, phi));
d2EEdTheta2    = simplify(Jacobian(dEEdTheta, Theta));
d2EEdPhi2      = simplify(Jacobian(dEEdPhi, phi));
d2EEdPhidTheta = simplify(Jacobian(dEEdPhi, Theta));
d2EEdThetadPhi = simplify(Jacobian(dEEdTheta, phi));

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
assume(V,'real');
for i=1:length(phi)
    for j=1:i
        sysVar = [sysVar, V(i,j)];
    end
end
% fundamental matrix
Phi = sym('Phi', [numel(phi), numel(phi)]);
assume(Phi,'real');
sysVar = [sysVar reshape(Phi,1,[])];
dVdt = Afunc*V + V*Afunc.' + subs(Efunc*Efunc.');
RHS = [S*reactionFlux; dVdt(find(triu(ones(size(dVdt))))); reshape(Afunc*Phi,[],1)];

systemJacobian = simplify(Jacobian(RHS, sysVar));
systemJacobian_diag = diag(systemJacobian); % just the diagonal

%% initial sensitivities
disp('computing initial sensitivities')

if COMPUTE_Y0
    % dY/dTheta
    S0  = simplify(Jacobian(Y0, Theta));
    % d2Y/dTheta2
    S20 = simplify(Jacobian(S0, Theta));
else
    S0  = zeros(1,nvar*npar); % placeholder
    S20 = zeros(1,nvar*npar*npar); % placeholder
end

if COMPUTE_V0
    % dV/dTheta
    SV0  = simplify(Jacobian(V0, Theta));
    % d2V/dTheta2
    S2V0 = simplify(Jacobian(SV0, Theta));
else
    SV0  = zeros(1,nvar*(nvar+1)/2*npar); % placeholder
    S2V0 = zeros(1,nvar*(nvar+1)/2*npar*npar); % placeholder
end

%% generate the C library
disp('Generating C source')

% old path
olddir = cd([dirName '/matlab']);

% requires phi, t, Theta
objs1 = {'reactionFlux', 'J', 'dFdTheta', 'd2fdTheta2', 'Afunc', 'dAdTheta', ...
         'dAdPhi', 'd2AdPhi2', 'd2AdTheta2', 'd2AdThetadPhi', 'd2AdPhidTheta', ...
         'd2AdPhidTheta', 'Efunc', 'EEfunc', 'dEEdTheta', 'd2EEdTheta2', 'dEEdPhi', ...
         'd2EEdPhi2', 'd2EEdPhidTheta', 'd2EEdThetadPhi', 'systemJacobian_diag'};

% requires only Theta
objs2 = {'S0','S20','SV0','S2V0','Y0','V0'};

addpath(lnaMatlabDir)
for i = 1:length(objs1)
    fprintf('Generating %s source...\n', objs1{i})
    genCCode(eval(objs1{i}), objs1{i}, {phi, t, Theta}, {'VECTOR','SCALAR','VECTOR'});
end

for i = 1:length(objs2)
    fprintf('Generating %s source...\n', objs2{i})
    genCCode(sym(eval(objs2{i})), objs2{i}, {Theta}, {'VECTOR'});  % make sure the zero matrices are symbolic vars
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
tmp1=num2cell(S*F);
tmp2=num2cell(phi);
tmp = solve(tmp1{:},tmp2{:}); % solve for phi

if isempty(tmp)
    warning('Could not solve for steady state initial conditions.  Please supply the initial condition explicitly.')
    Y0 = [];
    return
end

if isstruct(tmp)
    tmp2 = struct2cell(tmp);
else
    tmp2 = {tmp};
end
Y0 = [tmp2{:}];

end

function [V0] = solveSS_var(A,E,F,S,phi,Theta,Y0) %#codegen
% solve for the steady state covariance matrix
% set the fluctuations to steady-state

% dVdt = AV + VA' + EE'
V = sym('V', [numel(phi), numel(phi)]);
assume(V, 'real');
V2 = V.';
V = subs(V2,V2(find(triu(ones(size(V2)),1))),V(find(triu(ones(size(V)),1)))); % symmetrize

dVdt = A*V + V*A.' + subs(E*E.');

t = sym('t','real');
% set the upper triangular components to zero flux!

nvar = numel(phi);

f1 = matlabFunction(dVdt, 'vars', {phi, t, Theta, V}); % function handle to dVdt
f2 = f1(Y0,t,Theta,V); % matrix of symbols
f3 = f2(find(triu(ones(size(f2)))));
% convert to cell arrays for backwards compatibility
tmp1=num2cell(f3);
tmp2=num2cell(V(find(triu(V))));
V0 = solve(tmp1{:},tmp2{:});
if isstruct(V0)
    tmp= struct2cell(V0);
else
    tmp={V0};
end
V0 =[tmp{:}].';

if isempty(V0)
    warning('Could not solve for steady state variance.  Setting initial variance to zero.')
    V0 = zeros(nvar*(nvar+1)/2,1);
end

end

function J=Jacobian(x,y)
    J = jacobian(reshape(x,[],1),y);
end
