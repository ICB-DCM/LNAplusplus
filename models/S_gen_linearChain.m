function S = S_gen_linearChain(N)
% create a stoichiometric matrix for a linear chain of reactions
% with N species
% N-1 forward reactions and N-1 backward reaction

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
    