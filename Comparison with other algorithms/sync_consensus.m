%This function implements synchronous continuous consensus dynamics. The
%graph has adjacency matrix A and two communities of sizes V1 and V2. The
%community of a fraction eta of the nodes is known. The number of
%iterations is passed as an input.

function [duration, error] = sync_consensus(A, V1, V2, eta, iterations)
tic
%ALGORITHM INITIALIZATION
votes = zeros(V1 + V2, 1);
V1R = round(eta * V1);
V2R = round(eta * V2);
votes(1 : V1R) = ones(V1R, 1);
votes(V1 + 1 : V1 + V2R) = -ones(V2R, 1);
votes_0 = votes;
%Update matrices
A = A + A';
d = sum(A, 2);
D = diag([zeros(V1R, 1); 1 ./ d(V1R + 1 : V1); zeros(V2R, 1); 1 ./ d(V1 + V2R + 1 : V1 + V2)]);
A = D * A;
b = diag([ones(1, V1R) zeros(1, V1 - V1R) ones(1, V2R) zeros(1, V2 - V2R)]) * votes_0;
%ALGORITHM ITERATIONS
for k = 1 : iterations
    %Updating votes
    votes = A * votes + b;
end
error = (sum(votes(1 : V1) <= 0) + sum(votes(V1 + 1 : V1 + V2) >= 0)) / (V1 + V2);
duration = toc;
end