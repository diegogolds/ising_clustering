clear all

%PROBLEM PARAMETERS
a = 3;
b = 1;
v1 = 1;
v2 = 1;
n = 50000;
V1 = round(v1 * n);
V2 = round(v2 * n);
lambda = @(x) log(x);
p = a * lambda(n) / n;
q = b * lambda(n) / n;

%SIMULATION PARAMETERS
alpha = 10;
na = [v1 v2] / (v1 + v2);
cab = (V1 + V2) * [p q; q p];
bp_max_iter = 20;
is_max_iter = 20 * (V1 + V2) * lambda(n);

%OPTIONS
rounds = 10;
eta = 0.06;

%EXPERIMENTS
target_error = (0.4 : 0.4 : 2) / 100;
%target_error = (0.2 : 0.2 : 2) / 100;
edges = zeros(rounds, 1);
%Data vectors for our algorithm
is_iter = zeros(length(target_error), rounds);
is_erro = zeros(length(target_error), rounds);
is_flip = zeros(length(target_error), rounds);
is_fail = zeros(length(target_error), 1);
%Data vectors for Belief Propagation
bp_iter = zeros(length(target_error), rounds);
bp_erro = zeros(length(target_error), rounds);
bp_fail = zeros(length(target_error), 1);
for r = 1 : rounds
    %Stochastic block model
    tic
    disp("Stochastic block model")
    [A, ~, ~] = sbm([V1; V2], [p q; q p]);
    edges(r) = nnz(A);
    A = A + A';
    toc
    for k = 1 : length(target_error)
        disp("Round = " + r + " - Target error = " + target_error(k) * 100 + "%");
        %Our algorithm
        disp("Ising:");
        tic
        [error, flips, iterations] = ising(A, V1, V2, n, alpha, eta, lambda, target_error(k), is_max_iter);
        toc
        is_iter(k, r) = iterations;
        is_erro(k, r) = error;
        is_flip(k, r) = flips;
        if error > target_error(k)
            is_fail(k) = is_fail(k) + 1;
        end
        disp("- Error = " + 100 * error + "%");
        disp("- Normalized iterations = " + iterations / edges(r));
        %Semi-Supervised Belief Propagation using exact SBM parameters
        disp("Belief Propagation:");
        tic
        [error, iterations] = semi_supervised_bp(A, na, cab, eta, target_error(k), bp_max_iter);
        toc
        bp_iter(k, r) = iterations;
        bp_erro(k, r) = error;
        if error > target_error(k)
            bp_fail(k) = bp_fail(k) + 1;
        end
        disp("- Error = " + 100 * error + "%");
        disp("- Iterations = " + iterations);
    end
end

save("10_5_nodes_eta_06");