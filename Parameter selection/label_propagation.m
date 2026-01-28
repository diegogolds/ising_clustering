%This function implements the continuous-time version of label propagation.
%Most of the names of the variables follow the notation in the paper
%closely. The adjacency matrix of the graph is denoted by A. The maximum
%running time and tolerance with respect to the relative error as passed as
%arguments. The algorithm stops when the number of iterations exceeds the
%maximum value or the relative error is below the tolerance.

function [duration, error, iterations] = label_propagation(A, V1, V2, eta, max_iter, tolerance)
tic
%ALGORITHM INITIALIZATION
I = 1 : V1 + V2;
%Initial configuration
spins = sign(rand(1, V1 + V2) - 0.5);
for i = 1 : round(eta * V1)
    spins(i) = 1;
end
for i = V1 + 1 : V1 + round(eta * V2)
    spins(i) = -1;
end
%Energy variation
dE = 2 * spins .* (spins * (A + A'));
%Initial condition
k = 1;
X = sum(spins(1 : V1));
Y = sum(spins(V1 + 1 : V1 + V2));
%ALGORITHM ITERATIONS
error = (V1 - X + V2 + Y) / (2 * (V1 + V2));
J = I(dE < 0);
while ~isempty(J) && k < max_iter && error > tolerance
    %Selecting node that flips
    i = J(randi(length(J)));
    %Updating spins and energies
    spins(i) = -spins(i);
    aux = dE + 4 * spins(i) * spins .* [A(1 : i, i)' A(i, i + 1 : V1 + V2)];
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
    J = I(dE < 0);
end
iterations = k - 1;
duration = toc;
end

