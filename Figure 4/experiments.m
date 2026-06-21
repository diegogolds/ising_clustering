clear all

%PROBLEM PARAMETERS
a = 3;
b = 1;
v1 = 1;
v2 = 1;
n = 500000;
V1 = round(v1 * n);
V2 = round(v2 * n);
lambda = @(x) log(x);
p = a * lambda(n) / n;
q = b * lambda(n) / n;

%SIMULATION PARAMETERS
alpha = 10;
na = [v1 v2] / (v1 + v2);
cab = (V1 + V2) * [p q; q p];

%OPTIONS
eta = 0.02;
bp_max_iter = 6;
bp_iter = 1 : bp_max_iter;
is_step = 1000000;
is_max_step = 10;
is_iter = is_step : is_step : is_max_step * is_step;

%EXPERIMENTS
is_erro = zeros(1, is_max_step);
is_flip = zeros(1, is_max_step);
is_time = zeros(1, is_max_step);
%Data vectors for Belief Propagation
bp_erro = zeros(1, bp_max_iter);
bp_time = zeros(1, bp_max_iter);
%Stochastic block model
tic
disp("Stochastic block model")
[A, ~, ~] = sbm([V1; V2], [p q; q p]);
A = A + A';
edges = nnz(A);
toc
for k = 1 : bp_max_iter
    disp("Round = " + k);
    %Semi-Supervised Belief Propagation using exact SBM parameters
    disp("Belief Propagation:");
    tic
    [error, ~] = semi_supervised_bp(A, na, cab, eta, -1, bp_iter(k));
    time = toc;
    bp_erro(k) = error;
    bp_time(k) = time;
    disp("- Error = " + 100 * error + "%");
    disp("- Time = " + time);
end

for k = 1 : is_max_step
    disp("Round = " + k);
    %Our algorithm
    disp("Ising:");
    tic
    [error, flips, ~] = ising(A, V1, V2, n, alpha, eta, lambda, -1, k * is_step);
    time = toc;
    is_erro(k) = error;
    is_flip(k) = flips;
    is_time(k) = time;
    disp("- Error = " + 100 * error + "%");
    disp("- Time = " + time);
end

%SCORES
is_score = is_iter - is_flip + (3 + (a + b) * lambda(n)) * is_flip;
bp_score = (2 * (1 - eta)^2 * (a + b) * n * lambda(n) + 2 * (1 - eta) * n) * bp_iter;

save("data", "a", "b", "v1", "v2", "n", "lambda", "alpha", "eta", "is_iter", "is_erro", "is_flip", "is_time", "is_score", "bp_iter", "bp_erro", "bp_time", "bp_score");