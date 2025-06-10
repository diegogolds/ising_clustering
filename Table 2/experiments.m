clear all

%FIXED PARAMETERS
lambda = @(x) log(x);
n = 5000;
v1 = 1;
v2 = 1;
V1 = round(v1 * n);
V2 = round(v2 * n);
a = 3;
b = 1;

%OPTIONS
iter_per_node = 20;
rounds = 10;

%EXPERIMENTS
iterations1 = iter_per_node * (V1 + V2);
iterations2 = iter_per_node;
iterations3 = iter_per_node * (V1 + V2) * (a * lambda(n) + b * lambda(n)) / 2;
d = zeros(10, 8, rounds); %Real duration in seconds
e = zeros(10, 8, rounds); %Relative error
i = zeros(10, 8, rounds); %Number of iterations
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
        %alpha = 10 and beta = infty
        [duration, error, iterations] = ising_inf_beta(A, V1, V2, n, 10, eta, lambda, iterations1);
        d(k, 1, r) = duration;
        e(k, 1, r) = error;
        i(k, 1, r) = iterations;
        disp("-alpha = 10 and beta = infty: " + 100 * error + "% error in " + iterations + " iterations and " + duration + " sec");
        %Asynchronous Consensus
        [duration, error] = async_consensus(A, V1, V2, eta, iterations1);
        d(k, 2, r) = duration;
        e(k, 2, r) = error;
        i(k, 2, r) = iterations1;
        disp("-Asynchronous Consensus: " + 100 * error + "% error in " + iterations1 + " iterations and " + duration + " sec");
        %Synchronous Consensus
        [duration, error] = sync_consensus(A, V1, V2, eta, iterations2);
        d(k, 3, r) = duration;
        e(k, 3, r) = error;
        i(k, 3, r) = iterations2;
        disp("-Synchronous Consensus: " + 100 * error + "% error in " + iterations2 + " iterations and " + duration + " sec");
        %Gossiping
        [duration, error] = gossiping(I, J, V1, V2, eta, iterations3);
        d(k, 4, r) = duration;
        e(k, 4, r) = error;
        i(k, 4, r) = iterations3;
        disp("-Gossiping: " + 100 * error + "% error in " + iterations3 + " iterations and " + duration + " sec");
        %Standard Laplacian
        [duration, error] = gen_laplacian(0.95, 1, A, V1, V2, eta, iterations2);
        d(k, 5, r) = duration;
        e(k, 5, r) = error;
        i(k, 5, r) = iterations2;
        disp("-Standard Laplacian: " + 100 * error + "% error in " + iterations2 + " iterations and " + duration + " sec");
        %Normalized Laplacian
        [duration, error] = gen_laplacian(0.95, 0.5, A, V1, V2, eta, iterations2);
        d(k, 6, r) = duration;
        e(k, 6, r) = error;
        i(k, 6, r) = iterations2;
        disp("-Normalized Laplacian: " + 100 * error + "% error in " + iterations2 + " iterations and " + duration + " sec");
        %Page Rank
        [duration, error] = gen_laplacian(0.95, 0, A, V1, V2, eta, iterations2);
        d(k, 7, r) = duration;
        e(k, 7, r) = error;
        i(k, 7, r) = iterations2;
        disp("-Page Rank: " + 100 * error + "% error in " + iterations2 + " iterations and " + duration + " sec");
        %Poisson Learning
        [duration, error] = poisson(A, V1, V2, eta, iterations2);
        d(k, 8, r) = duration;
        e(k, 8, r) = error;
        i(k, 8, r) = iterations2;
        disp("-Poisson Learning: " + 100 * error + "% error in " + iterations2 + " iterations and " + duration + " sec");
    end
    save("k_" + k);
end

eta = 0.01 : 0.01 : 0.1;
algorithms = ['\alpha = 10 and \beta = \infty', 'Asynchronous Consensus', 'Synchronous Consensus', 'Gossiping', 'Standard Laplacian', 'Normalized Laplacian', 'Page Rank', 'Poisson Learning'];
save("log_3");