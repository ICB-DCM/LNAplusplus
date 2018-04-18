function [ S, Sed, P ] = SBML2StoichProp( sbmlModel, doReplacePower )

    SpeciesIDs = {sbmlModel.species.id};
    Ss = length(sbmlModel.species);
    Rs = length(sbmlModel.reaction);
    
    S = zeros(Ss, Rs);
    Sed = zeros(Ss, Rs);
    P = cell(Rs,1);
    for rIdx = 1:Rs
        reaction = sbmlModel.reaction(rIdx);
        [found,eductLoc]=ismember({reaction.reactant.species}, SpeciesIDs);
        if (any(~found))
            error('Unknown reactant in reaction %d', rIdx);
        end
        if (~isempty(eductLoc))
            S(eductLoc,rIdx) = S(eductLoc,rIdx)-[reaction.reactant.stoichiometry]';
            Sed(eductLoc,rIdx) = Sed(eductLoc,rIdx)-[reaction.reactant.stoichiometry]';
        end
        
        [found,productLoc]=ismember({reaction.product.species}, SpeciesIDs);
        if (any(~found))
            error('Unknown product in reaction %d', rIdx);
        end
        if (~isempty(productLoc))
            S(productLoc,rIdx) = S(productLoc,rIdx) + [reaction.product.stoichiometry]';
        end
       
        P{rIdx} = reaction.kineticLaw.formula;
        
        if (doReplacePower)
            P{rIdx} = regexprep(P{rIdx}, 'power\(([\w_]+),(\d+))', '$1^$2');
        end
        
%         if (doAdaptHomodimerProp)
%             P{rIdx} = regexprep(P{rIdx}, '([\w_]+)\^2([^\d]|$)', '(($1)*($1-1)/2)');
%             P{rIdx} = regexprep(P{rIdx}, '([\w_]+)\^3([^\d]|$)', '(($1)*($1-1)*($1-2)/6)');
%             if (~isempty(regexp(P{rIdx}, '[\w_]+\^\d','ONCE')))
%                 error('Cannot adapt propensity for 3+ order reaction');
%             end
%         end        
    end
    
    D = zeros(Rs,Rs);
    for r = 1:Rs % for each reaction
        for s = 1:Ss % for each species
            D(r,:) = D(r,:) | (Sed(s,:)&S(s,r));
        end
    end
%     ReacDep = cell(Rs,1);
%     for r = 1:Rs
%         ReacDep{r} = find(D(r,:));
%     end    

    % replace all species names with phi(index) and parameter names with
    % Theta(index) in the propensity function expressions
    
    Ts = length(sbmlModel.parameter);
    P = cell2sym(P);
    
    subsPhi = arrayfun(@(k)sprintf('phi(%d)',k),1:Ss,'UniformOutput',false);
    subsTheta = arrayfun(@(k)sprintf('Theta(%d)',k),1:Ts,'UniformOutput',false);

    P = subs(P, {sbmlModel.species.name}, subsPhi);
    P = subs(P, {sbmlModel.parameter.name}, subsTheta);
    
end


