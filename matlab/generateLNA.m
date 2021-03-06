% Generate the executable and wrapper function for the specified SBML model
% Input:
% sbmlFileName: the path of the SBML model specification file
% modelName: the desired name of the model 
% computeSS:  this argument is optional and can have only one of the values ?Y0?, ?V0?, ?BOTH?, or ?NONE?.
%  Y0:  LNA++ will attempt to calculate the steady state value of the MRE of the stochastic system specified.  If Y0 is not explicitly provided when invoking the model, the computed value will be automatically substituted. 
%  V0: LNA++ will attempt to calculate the steady state value of the variance of the stochastic system specified.  If V0 is not explicitly provided when invoking the model, the computed value will be automatically substituted.
%  BOTH: LNA++ will compute both V0 and Y0 if possible.  If Y0 or V0 is not explicitly provided when invoking the model, the computed values will be automatically substituted.
%  NONE: LNA++ will not compute any steady state values.  Invocation of the model without explicitly specifying the initial conditions for Y0 and V0 will result in an error message being generated.
% varargin: additional arguments passed to compileLNA, see there
function generateLNA(sbmlFileName, modelName, computeSS, varargin)

sbmlModel = TranslateSBML(sbmlFileName);
[ S, ~, P ] = SBML2StoichProp(sbmlModel, true);

orig_state = warning('off', 'symbolic:generate:FunctionNotVerifiedToBeValid');
f = matlabFunction(P, 'Vars', {'phi','t','Theta'});
warning(orig_state)

% species
phi = cell2sym({sbmlModel.species.name});
Theta = cell2sym({sbmlModel.parameter.name});
    
%% flux vector

% generate the C code for the components of the LNA model
generateLNAComponents(modelName, S, f, phi, Theta, computeSS);
% compile code
if ~isempty(varargin)
    compileLNA(modelName, varargin); 
else
    compileLNA(modelName)
end
