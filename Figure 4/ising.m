%This function implements Algorithm 1 of https://arxiv.org/pdf/2506.09223
%when the inverse temperature beta is infinite.
%A            = adjacency matrix (sparse)
%V1           = number of nodes in community 1
%V2           = number of nodes in community 2
%n            = scaling parameter
%alpha        = quadratic penalty
%eta          = fraction of revealed community labels
%lambda       = average degree
%target_error = relative error for stopping
%max_iter     = maximum number of iterations

function [error, flips, iterations] = ising(A, V1, V2, n, alpha, eta, lambda, target_error, max_iter)
%ALGORITHM INITIALIZATION
alpha_n = alpha * lambda(n) / n;
%Uniformly random initial spins
spin = 2 * (rand(1, V1 + V2) < 0.5) - 1;
%Revealed community labels
spin(1 : round(eta * V1)) = ones(1, round(eta * V1));
spin(V1 + 1 : V1 + round(eta * V2)) = -ones(1, round(eta * V2));
%Magnetization over each community
X = sum(spin(1 : V1));
Y = sum(spin(V1 + 1 : V1 + V2));
%Magnetization over each neighborhood
M = spin * A;
%ALGORITHM ITERATIONS
error = (V1 - X + V2 + Y) / (2 * (V1 + V2));
flips = 0;
k = 0;
while error > target_error && k < max_iter
    %Selecting node uniformly at random
    i = randi(V1 + V2);
    %Evaluating sign of energy variation
    dE = spin(i) * (M(i) - alpha_n * (X + Y)) + alpha_n;
    if dE < 0
        flips = flips + 1;
        %Updating spin
        spin(i) = -spin(i);
        %Updating neighborhood magnetizations
        neighbors = find(A(:, i));
        M(neighbors) = M(neighbors) + 2 * spin(i);
        %Updating community magnetizations
        if i <= V1
            X = X + 2 * spin(i);
        else
            Y = Y + 2 * spin(i);
        end
    end
    %Further updates
    error = (V1 - X + V2 + Y) / (2 * (V1 + V2));
    k = k + 1;
end
iterations = k - 1;
end

