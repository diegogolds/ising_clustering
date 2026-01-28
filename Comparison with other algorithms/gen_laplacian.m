%This function implements generalized Laplacian learning as in the paper
%by Avrachenkov, Goncalves, Mishenin and Sokol. Here A is the adjacency
%matrix of the graph (upper triangular and sparse), V1 and V2 are the
%community sizes, eta is the fraction of revealed labels per community and
%iterations is the number of iterations.

function [duration, error] = gen_laplacian(alpha, sigma, A, V1, V2, eta, iterations)
tic
%INITIALIZATION
A = A + A';
deg = sum(A, 2);
D1 = diag(deg .^ (-sigma));
D2 = diag(deg .^ (sigma - 1));
%Revealed labels
R1 = round(eta * V1);
R2 = round(eta * V2); 
Y = [ones(R1, 1) zeros(R1, 1); zeros(V1 + V2 - R1 - R2, 2); zeros(R2, 1) ones(R2, 1)];
%ALGORITHM ITERATIONS;
U = zeros(V1 + V2, 2);
for k = 1 : iterations
    U = alpha * D1 * A * D2 * U + (1 - alpha) * Y;
end
error = (sum(U(1 : V1, 2) > U(1 : V1, 1)) + sum(U(V1 + 1 : V1 + V2, 1) > U(V1 + 1 : V1 + V2, 2))) / (V1 + V2);
duration = toc;
end