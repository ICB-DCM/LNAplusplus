function varargout=logLikelihood(z, options)

z=exp(z);

g_m=z(1);
g_p=z(2);

data = options.data;

k_m=20;
k_p=25;
merr=0;
y0=[2 200];

P = 0;
grad = zeros(1,length(z)+1);

for n=1:1 %length(data)
    data_n = data{n}(:,3);
    tspan = data{n}(:,1);
    if nargout==1
        [mu, var] = BirthDeath_LNA([k_m, k_p, g_m, g_p], tspan, 2, merr, y0);

    else
%         [mu,var,SensMRE,SensVar] = BirthDeath_LNA(tspan, [k_m, k_p, g_m, g_p], merr, 2, y0);
        [mu,var,SensMRE,SensVar] = BirthDeath_LNA([k_m, k_p, g_m, g_p], tspan, 2, merr, y0);

        % todo: add derivative of likelihood w.r.t. g_m and g_p
        SigmaInv    = eye(size(var))/var;
        g = gradient_ll(data{n}, mu, z, SensVar, SigmaInv, SensMRE) .* [z' merr];
        grad = grad + g;
    end
    
    P = P + sum(logmvnpdf(data_n', mu, var));
end

if strcmp(options.sign,'negative')
    P=-P;
end

if nargout>0
    varargout{1}=P;
end

if nargout>1
%     grad           = grad(PLEoptions.grad_ind);   
    if strcmp(options.sign,'negative')
        grad=-grad;
    end
    varargout{2} = grad(1:length(z));
end
    
end


function g = gradient_ll(data, Y, par, dSigma, SigmaInv, Sens)
% return the grad of the negative log-likelihood

% mu = Y(2,:)';
mu=Y;

% x       = data.int'*lambda;
x       = data(:,3)';
npar    = length(par);
g       = zeros(1,npar+1);
dSigma  = squeeze(dSigma);
ixPar   = [3 4];

for ix=1:npar
    k=ixPar(ix);
    dmu_dk = squeeze(Sens(k,:)); % compute in linear space
    g(ix) = -0.5*trace(SigmaInv*dSigma(:,:,k)) + ...
        0.5*(dmu_dk)*SigmaInv*(x-mu)' + 0.5*(x-mu)*(SigmaInv*dSigma(:,:,k)*SigmaInv)*(x-mu)' + 0.5*(x-mu)*SigmaInv*(dmu_dk');
end

%df / dmerr 
g(npar+1) = -0.5*trace(SigmaInv) + 0.5*(x-mu)*SigmaInv*SigmaInv*(x-mu)';

end