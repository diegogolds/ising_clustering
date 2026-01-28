%This function implements asynchronous continuous consensus dynamics. The
%graph has adjacency matrix A and two communities of sizes V1 and V2. The
%community of a fraction eta of the nodes is known. The number of
%iterations is passed as input

function [duration, error] = async_consensus(A, V1, V2, eta, iterations)
tic
%ALGORITHM INITIALIZATION
votes = zeros(V1 + V2, 1);
V1R = round(eta * V1);
V2R = round(eta * V2);
votes(1 : V1R) = ones(V1R, 1);
votes(V1 + 1 : V1 + V2R) = -ones(V2R, 1);
%Update matrices
A = A + A';
d = sum(A, 2);
%ALGORITHM ITERATIONS
for k = 1 : iterations
    %Selecting node uniformly at random
    if rand < (V1 - V1R) / (V1 + V2 - V1R - V2R)
        i = V1R + randi(V1 - V1R);
    else
        i = V1 + V2R + randi(V2 - V2R);
    end
    %Updating vote
    votes(i) = (votes(i) + A(i, :) * votes) / (d(i) + 1);
end
error = (sum(votes(1 : V1) < 0) + sum(votes(V1 + 1 : V1 + V2) > 0)) / (V1 + V2);
duration = toc;
end