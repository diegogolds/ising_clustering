%This function implements the continuous-time version of the community
%detection algorithm based on the Glauber dynamics for the Ising model.
%Most of the names of the variables follow the notation in the paper
%closely. The adjacency matrix of the graph is denoted by A. The maximum
%running time and tolerance with respect to the relative error as passed as
%arguments. The algorithm stops when the running time exceeds the maximum
%running time or the relative error is smaller than the tolerance.

function [T, X, Y] = ising_fin_beta(A, V1, V2, n, alpha, beta, eta, lambda, time, tolerance, scaling, progressbar)
%ALGORITHM INITIALIZATION
tic
disp("Initializing the algorithm")
%Parameters
if scaling == 0
    alpha_n = alpha * lambda(n) / n;
    beta_n = beta / lambda(n);
else
    alpha_n = alpha * lambda(n) / n;
    beta_n = beta;
end
%Initial configuration
spins = sign(rand(1, V1 + V2) - 0.5);
for i = 1 : V1
    if rand < abs(eta(1))
        spins(i) = sign(eta(1));
    end
end
for i = V1 + 1 : V1 + V2
    if rand < abs(eta(2))
        spins(i) = sign(eta(2));
    end
end
%Energy variation
c = ones(1, V1 + V2);
c(1 : V1) = n / V1;
c(V1 + 1 : V1 + V2) = n / V2;
cspins = c .* spins;
dE = 2 * spins .* (spins * (A + A')) - 2 * spins * alpha_n * sum(cspins) + 2 * alpha_n * c;
%Data vectors
events = time * (V1 + V2);
T = zeros(events, 1);
X = zeros(events, 1);
Y = zeros(events, 1);
%Initial condition
k = 1;
t = 0;
X(k) = sum(spins(1 : V1));
Y(k) = sum(spins(V1 + 1 : V1 + V2));
%Initializing waitbar
if progressbar == 1
    progress = 0;
    w = waitbar(progress, "Progress");
end
toc
%ALGORITHM ITERATIONS
tic
disp("Running the algorithm")
%Dynamics
error = (V1 - X(k) + V2 + Y(k)) / (2 * (V1 + V2));
while t < time && error > tolerance
    %Computing time until next flip
    [dt, i] = min(exprnd(1 + exp(beta_n * dE)));
    %Updating spins and energies
    spins(i) = -spins(i);
    cspins(i) = -cspins(i);
    aux = dE + 4 * spins(i) * spins .* [A(1 : i, i)' A(i, i + 1 : V1 + V2)] - 4 * alpha_n * cspins(i) * spins;
    aux(i) = -dE(i);
    dE = aux;
    %Updating data vectors
    k = k + 1;
    t = t + dt;
    T(k) = t;
    if i < V1 + 1
        if spins(i) > 0
            X(k) = X(k - 1) + 2;
        else
            X(k) = X(k - 1) - 2;
        end
        Y(k) = Y(k - 1);
    else
        if spins(i) > 0
            Y(k) = Y(k - 1) + 2;
        else
            Y(k) = Y(k - 1) - 2;
        end
        X(k) = X(k - 1);
    end
    %Updating error
    error = (V1 - X(k) + V2 + Y(k)) / (2 * (V1 + V2));
    %Updating waitbar
    if progressbar == 1
        if t / time - progress > 0.1
            progress = progress + 0.1;
            waitbar(progress, w, "Progress");
        end
    end
end
%Closing waitbar
if progressbar == 1
    close(w);
end
T = [T(1 : k - 1); time];
X = [X(1 : k - 1); X(k - 1)];
Y = [Y(1 : k - 1); Y(k - 1)];
toc
end

