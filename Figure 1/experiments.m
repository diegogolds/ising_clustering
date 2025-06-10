clear all
close all

%GRAPH PARAMETERS
lambda = @(x) sqrt(log(x));
n = 5000;
v1 = 1;
v2 = 1;
a = 5;
b = 1;

%ALGORITHM PARAMETERS
eta = 0.1;
alpha = 10;
beta = 1;

%SIMULATION PARAMETERS
progressbar = 0;
tolerance = -1;
rounds = 100;
scaling = 1;
time = 5;

%DATA VECTORS
t = 0 : time / 500 : time;
Xfin = zeros(rounds, length(t));
Yfin = zeros(rounds, length(t));
Xinf = zeros(rounds, length(t));
Yinf = zeros(rounds, length(t));
Xzrf = zeros(rounds, length(t));
Yzrf = zeros(rounds, length(t));
Xzri = zeros(rounds, length(t));
Yzri = zeros(rounds, length(t));

%SIMULATIONS
V1 = round(v1 * n);
V2 = round(v2 * n);
p = a * lambda(n) / n;
q = b * lambda(n) / n;
for r = 1 : rounds
    disp("Round " + r)
    %STOCHASTIC BLOCK MODEL
    tic
    disp("Stochastic block model")
    A = sbm([V1; V2], [p q; q p]);
    toc
    %FINITE INVERSE TEMPERATURE
    disp('Fininte inverse temperature')
    [T1, X1, Y1] = ising_fin_beta(A, V1, V2, n, alpha, beta, [eta -eta], lambda, time, tolerance, scaling, progressbar);
    %INFINITE INVERSE TEMPERATURE
    disp('Infinite inverse temperature')
    [T2, X2, Y2] = ising_inf_beta(A, V1, V2, n, alpha, [eta -eta], lambda, time, tolerance, progressbar);
    %FINITE INVERSE TEMPERATURE WITHOUT ALPHA
    disp('Fininte inverse temperature')
    [T3, X3, Y3] = ising_fin_beta(A, V1, V2, n, 0, beta, [eta -eta], lambda, time, tolerance, scaling, progressbar);
    %INFINITE INVERSE TEMPERATURE WITHOUT ALPHA
    disp('Fininte inverse temperature')
    [T4, X4, Y4] = ising_inf_beta(A, V1, V2, n, 0, [eta -eta], lambda, time, tolerance, progressbar);
    tic
    disp('Sampling')
    for i = 1 : length(t)
        j1 = find(T1 >= t(i), 1);
        Xfin(r, i) = X1(j1);
        Yfin(r, i) = Y1(j1);
        j2 = find(T2 >= t(i), 1);
        Xinf(r, i) = X2(j2);
        Yinf(r, i) = Y2(j2);
        j3 = find(T3 >= t(i), 1);
        Xzrf(r, i) = X3(j3);
        Yzrf(r, i) = Y3(j3);
        j4 = find(T4 >= t(i), 1);
        Xzri(r, i) = X4(j4);
        Yzri(r, i) = Y4(j4);
    end
    toc
end
save("n5000a5eta01sqrtlog");