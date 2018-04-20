% S_gen_linearChain.m creates the stoichiometric matrix for a linear chain 
% of reactionswith N species

function S = S_gen_linearChain(N)

S = zeros(N, 2*N);

% birth/death
S(1,1:2) = [1,-1];

% forward reactions
for i=1:N-1
    S(i:i+1,2+i) = [-1;1];
end

% backward reactions
for i=1:N-1
    S(i:i+1,N+1+i) = [1;-1];
end
    