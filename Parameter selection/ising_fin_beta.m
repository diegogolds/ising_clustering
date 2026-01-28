%This function implements the continuous-time version of the community
%detection algorithm based on the Glauber dynamics for the Ising model.
%Most of the names of the variables follow the notation in the paper
%closely. The adjacency matrix of the graph is denoted by A. The maximum
%running time and tolerance with respect to the relative error as passed as
%arguments. The algorithm stops when the number of iterations exceeds the
%maximum value or the relative error is smaller than the tolerance.

function [duration, error, iterations] = ising_fin_beta(A, V1, V2, n, alpha, beta, eta, lambda, max_iter, tolerance)
tic
%ALGORITHM INITIALIZATION
alpha_n = alpha * lambda(n) / n;
beta_n = beta;
%Initial configuration
spins = sign(rand(1, V1 + V2) - 0.5);
for i = 1 : round(eta * V1)
    spins(i) = 1;
end
for i = V1 + 1 : V1 + round(eta * V2)
    spins(i) = -1;
end
%Energy variation
dE = 2 * spins .* (spins * (A + A')) - 2 * spins * alpha_n * sum(spins) + 2 * alpha_n;
%Initial condition
k = 1;
X = sum(spins(1 : V1));
Y = sum(spins(V1 + 1 : V1 + V2));
%ALGORITHM ITERATIONS
error = (V1 - X + V2 + Y) / (2 * (V1 + V2));
while k < max_iter && error > tolerance
    %Computing time until next flip
    [~, i] = min(exprnd(1 + exp(beta_n * dE)));
    %Updating spins and energies
    spins(i) = -spins(i);
    aux = dE + 4 * spins(i) * spins .* [A(1 : i, i)' A(i, i + 1 : V1 + V2)] - 4 * alpha_n * spins(i) * spins;
    aux(i) = -dE(i);
    dE = aux;
    %Updating data vectors
    k = k + 1;
    if i < V1 + 1
        if spins(i) > 0
            X = X + 2;
        else
            X = X - 2;
        end
    else
        if spins(i) > 0
            Y = Y + 2;
        else
            Y = Y - 2;
        end
    end
    %Updating error
    error = (V1 - X + V2 + Y) / (2 * (V1 + V2));
end
iterations = k - 1;
duration = toc;
end

