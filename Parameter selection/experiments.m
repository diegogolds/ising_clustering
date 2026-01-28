clear all

%FIXED PARAMETERS
lambda = @(x) log(x);
n = 10000;
v1 = 1;
v2 = 0.75;
V1 = round(v1 * n);
V2 = round(v2 * n);
a = 10;
b = 1;

%OPTIONS
tolerance = -1;
max_iter = 5 * n;
rounds = 10;

%EXPERIMENTS
d = zeros(10, 4, rounds); %Real duration in seconds
e = zeros(10, 4, rounds); %Relative error
i = zeros(10, 4, rounds); %Number of iterations
for k = 1 : 10
    eta = 0.01 * k;
    for r = 1 : rounds
        disp("eta = " + eta + ", round = " + r);
        %Stochastic block model
        tic
        disp("Stochastic block model")
        p = a * lambda(n) / n;
        q = b * lambda(n) / n;
        [A, I, J] = sbm([V1; V2], [p q; q p]);
        toc
        %alpha = 0 and beta = 1
        [duration, error, iterations] = ising_fin_beta(A, V1, V2, n, 0, 1, eta, lambda, max_iter, tolerance);
        d(k, 1, r) = duration;
        e(k, 1, r) = error;
        i(k, 1, r) = iterations;
        disp("-alpha = 0 and beta = 1: " + 100 * error + "% error in " + iterations + " iterations and " + duration + " sec");
        %alpha = 0 and beta = infty
        [duration, error, iterations] = label_propagation(A, V1, V2, eta, max_iter, tolerance);
        d(k, 2, r) = duration;
        e(k, 2, r) = error;
        i(k, 2, r) = iterations;
        disp("-alpha = 0 and beta = infty: " + 100 * error + "% error in " + iterations + " iterations and " + duration + " sec");
        %alpha = 6 and beta = 1
        [duration, error, iterations] = ising_fin_beta(A, V1, V2, n, 6, 1, eta, lambda, max_iter, tolerance);
        d(k, 3, r) = duration;
        e(k, 3, r) = error;
        i(k, 3, r) = iterations;
        disp("-alpha = 6 and beta = 1: " + 100 * error + "% error in " + iterations + " iterations and " + duration + " sec");
        %alpha = 6 and beta = infty
        [duration, error, iterations] = ising_inf_beta(A, V1, V2, n, 6, eta, lambda, max_iter, tolerance);
        d(k, 4, r) = duration;
        e(k, 4, r) = error;
        i(k, 4, r) = iterations;
        disp("-alpha = 6 and beta = infty: " + 100 * error + "% error in " + iterations + " iterations and " + duration + " sec");
    end
    save("k_" + k);
end

eta = 0.01 : 0.01 : 0.1;
algorithms = ['\alpha = 0 and \beta = 1', '\alpha = 0 and \beta = \infty', '\alpha = 6 and \beta = 1', '\alpha = 6 and \beta = \infty'];
save("log_10");