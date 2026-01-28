%This function implements Poisson learning as in Algorithm 1 of the paper
%by Calder, Cook, Thorpe and Slepcev. Here A is the adjacency matrix of the
%graph (upper triangular and sparse), V1 and V2 are the community sizes,
%eta is the fraction of revealed labels per community and iterations is the
%number of iterations.

function [duration, error] = poisson(A, V1, V2, eta, iterations)
tic
%INITIALIZATION
W = A + A';
deg = sum(W, 2);
%Scalars
R1 = round(eta * V1);
R2 = round(eta * V2);
m = R1 + R2;
n = V1 + V2;
%Matrices
B = [(R2 / m) * ones(1, R1) zeros(1, n - R1); zeros(1, n - R2) (R1 / m) * ones(1, R2)];
Dinv = diag(deg .^ (-1));
D = diag(deg);
L = D - W;
%ALGORITHM ITERATIONS;
U = zeros(n, 2);
for k = 1 : iterations
    U = U + Dinv * (B' - L * U);
end
%Magnetization per community
error = (sum(U(1 : V1, 1) < U(1 : V1, 2)) + sum(U(V1 + 1 : V1 + V2, 1) > U(V1 + 1 : V1 + V2, 2))) / n;
duration = toc;
end