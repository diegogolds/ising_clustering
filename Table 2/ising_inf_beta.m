%This function implements the continuous-time version of label propagation.
%Most of the names of the variables follow the notation in the paper
%closely. The adjacency matrix of the graph is denoted by A. The number of
%iterations is passed as input.

function [duration, error, iterations] = ising_inf_beta(A, V1, V2, n, alpha, eta, lambda, max_iter)
tic
%ALGORITHM INITIALIZATION
I = 1 : V1 + V2;
alpha_n = alpha * lambda(n) / n;
%Initial configuration
spins = sign(rand(1, V1 + V2) - 0.5);
V1R = round(eta * V1);
V2R = round(eta * V2);
spins(1 : V1R) = ones(1, V1R);
spins(V1 + 1 : V1 + V2R) = -ones(1, V2R);
%Energy variation
dE = 2 * spins .* (spins * (A + A')) - 2 * spins * alpha_n * sum(spins) + 2 * alpha_n;
%ALGORITHM ITERATIONS
k = 0;
J = I(dE < 0);
while ~isempty(J) && k < max_iter
    %Computing time until next flip
    i = J(randi(length(J)));
    %Updating spins and energies
    spins(i) = -spins(i);
    aux = dE + 4 * spins(i) * spins .* [A(1 : i, i)' A(i, i + 1 : V1 + V2)] - 4 * alpha_n * spins(i) * spins;
    aux(i) = -dE(i);
    dE = aux;
    %Further updates
    k = k + 1;
    J = I(dE < 0);
end
error = (sum(spins(1 : V1) < 0) + sum(spins(V1 + 1 : V1 + V2) > 0)) / (V1 + V2);
iterations = k;
duration = toc;
end

