function varargout = BirthDeath_LNA(par, tspan, varargin)
if nargout > 4
    [Y Sigma Sens_MRE Sens_Var Sens2_MRE Sens2_Var] =  BirthDeath_LNA_mex(par, tspan, varargin);                   
    Y       = squeeze(Y);                                                                                  
    Sigma   = squeeze(Sigma);                                                                              
    Sens_MRE  = squeeze(Sens_MRE);                                                                         
    Sens_Var  = squeeze(Sens_Var);                                                                         
    Sens2_MRE = squeeze(Sens2_MRE);                                                                        
    Sens2_Var = squeeze(Sens2_Var);                                                                        
    varargout = {Y, Sigma, Sens_MRE, Sens_Var, Sens2_MRE, Sens2_Var};                                      
elseif nargout > 2                                                                                         
    [Y Sigma Sens_MRE Sens_Var] =  BirthDeath_LNA_mex(par, tspan, varargin);                                       
    Y       = squeeze(Y);                                                                                  
    Sigma   = squeeze(Sigma);                                                                              
    Sens_MRE  = squeeze(Sens_MRE);                                                                         
    Sens_Var  = squeeze(Sens_Var);                                                                         
    varargout = {Y, Sigma, Sens_MRE, Sens_Var};                                                            
else                                                                                                       
    [Y Sigma] =  BirthDeath_LNA_mex(par, tspan, varargin);                                                         
    Y       = squeeze(Y);                                                                                  
    Sigma   = squeeze(Sigma);                                                                              
    varargout = {Y, Sigma};                                                                                
end                                                                                                        
                                                                                                           
