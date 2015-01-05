% compileLNA.m
% (C) Justin Feigelman
% 2014
% Helmholtz Zentrum Muenchen

% convert all of the auto-generated functions into C code 
% and produce the Matlab executable
% input:
%  model: name of the model
%  S: the stoichiometric matrix, same as used for generating the m-files
%  Npar: the number of model parameters
% output executable is stored in ./mex directory

function compileLNA(model, S, Npar)

%% compile C code
disp('Compiling autogenerated code')
C_source = dir([model '/C/*.c']);
cmd_compile = [ sprintf(['mex -v CFLAGS=''$CFLAGS -Wno-header-guard -Wno-parentheses -Wno-incompatible-pointer-types'' -c -outdir %s/C ' ...
    '-I%s/C '], model, model), ...
    sprintf([model '/C/%s '], C_source.name)];

eval(cmd_compile)

%% compile CPP code
disp('Compiling LNA code')
% 
cmd_compile_cpp = ['mex CXXFLAGS=''$CXXFLAGS -Wno-header-guard -Wno-parentheses -Wno-incompatible-pointer-types'' -c ../src/computeLinearNoise.cpp ../src/mexLNA.cpp -I' model '/C -I/usr/local/include -I../include '];
%disp(cmd_compile_cpp)
eval( cmd_compile_cpp )

%% link executable
disp('Linking')
Cobj = dir([model '/C/*.o']);

% libstdcpp = ' /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib/libstdc++.dylib';
cmd_link = [ sprintf('mex -cxx -output %s/%s_LNA_mex  computeLinearNoise.o LNA_Onset.o %s', ...
    model, model, strjoin(fullfile(model, 'C', {Cobj.name}'))) ... libstdcpp ...
    ' -lsundials_cvodes -lsundials_nvecserial -lblitz -lstdc++ -lgsl -lgslcblas']; 
    
% disp(cmd_link);
eval(cmd_link)
% eval('delete *.o')
disp('Done')
%% generate wrapper
disp('Generating wrapper function')
% Npar = length(Theta);
% be supplied
wrapper_str = sprintf([
'%%Wrapper function auto-generated by LNA++\n'...
'%%Usage:\n'...
'%%\t[MRE,Var,dMREdTheta,dVardTheta,dMRE2dTheta2,dVar2dTheta2] = %s_LNA(Theta, tspan, [obsVar], [merr], [Y0], [V0])\n'...
'%%\tTheta: vector of model parameters\n'...
'%%\ttspan: vector of output times\n'...
'%% Optional arguments:\n'...
'%%\tobsVar: vector of indices (beginning with 1) corresponding to the species to output\n'...
'%%\tmerr: either a single number, or a vector the length of obsVar corresponding to the measurement error of each species\n'...
'%%\tY0: initial values for each species\n'...
'%%\tV0: initial values for each entry of the (upper triangular portion of the) covariance matrix\n'...
'%%\tin the default Matlab (i.e. column major) ordering\n'...
'%%Output:\n'...
'%%\tMRE: solution of the macroscopic rate equation for the observed variables specified\n' ...
'%%\tVar: block cross-covariance matrix for the observed variables specified\n'...
'%%\tdMREdTheta: sensitivity tensor of the MRE for each observed variable with respect to each parameter, at each time point\n' ...
'%%\tdVardTheta: sensitivity tensor of the block cross-covariance for each entry of the covariance matric w.r.t. all parameters\n'...
'%%\td2MREdTheta2: second order sensitivity tensor of the MRE solution\n'...
'%%\td2VardTheta2: second order sensitivity tensor of the block cross-covariance solution\n'...
'%%The algorithm uses the number of output arguments to decide whether or not to compute sensitivities, so it is less computationally\n'...
'%%intensive if fewer output arguments are specified\n'...
'function varargout = %s_LNA(par, tspan, varargin)\n' ...'addpath(''models/%s/matlab'')\n'...'addpath(''models/%s/mex'')\n'...
'if nargout > 4\n'...
'    [Y Sigma Sens_MRE Sens_Var Sens2_MRE Sens2_Var] =  %s_LNA_mex(par, tspan, varargin);                   \n'...
'    Y       = squeeze(Y);                                                                                  \n'...
'    Sigma   = squeeze(Sigma);                                                                              \n'...
'    Sens_MRE  = squeeze(Sens_MRE);                                                                         \n'...
'    Sens_Var  = squeeze(Sens_Var);                                                                         \n'...
'    Sens2_MRE = squeeze(Sens2_MRE);                                                                        \n'...
'    Sens2_Var = squeeze(Sens2_Var);                                                                        \n'...
'    varargout = {Y, Sigma, Sens_MRE, Sens_Var, Sens2_MRE, Sens2_Var};                                      \n'...
'elseif nargout > 2                                                                                         \n'...
'    [Y Sigma Sens_MRE Sens_Var] =  %s_LNA_mex(par, tspan, varargin);                                       \n'...
'    Y       = squeeze(Y);                                                                                  \n'...
'    Sigma   = squeeze(Sigma);                                                                              \n'...
'    Sens_MRE  = squeeze(Sens_MRE);                                                                         \n'...
'    Sens_Var  = squeeze(Sens_Var);                                                                         \n'...
'    varargout = {Y, Sigma, Sens_MRE, Sens_Var};                                                            \n'...
'else                                                                                                       \n'...
'    [Y Sigma] =  %s_LNA_mex(par, tspan, varargin);                                                         \n'...
'    Y       = squeeze(Y);                                                                                  \n'...
'    Sigma   = squeeze(Sigma);                                                                              \n'...
'    varargout = {Y, Sigma};                                                                                \n'...
'end                                                                                                        \n'...
'                                                                                                           \n'], ...
    model, model, model, model, model);

% write wrapper to a file
f = fopen(sprintf('%s/%s_LNA.m',model,model), 'w');

fprintf(f, '%s', wrapper_str);
fclose(f);
disp('Done')

end
