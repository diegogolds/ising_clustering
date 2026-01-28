%This function implements the Gossiping algorithm. The graph has adjacency
%matrix A and two communities of sizes V1 and V2. The community of a
%fraction eta of the nodes is known. The number of iterations is given.

function [duration, error] = gossiping(I, J, V1, V2, eta, iterations)
tic
%ALGORITHM INITIALIZATION
edges = length(I);
votes = zeros(V1 + V2, 1);
V1R = round(eta * V1);
V2R = round(eta * V2);
votes(V1 - V1R + 1 : V1) = ones(V1R, 1);
votes(V1 + 1 : V1 + V2R) = -ones(V2R, 1);
%ALGORITHM ITERATIONS
for k = 1 : iterations
    %Updating votes
    i = randi(edges);
    old1 = votes(I(i));
    old2 = votes(J(i));
    if (I(i) <= V1 - V1R) || (I(i) > V1 + V2R)
        new1 = (old1 + old2) / 2;
    else
        new1 = old1;
    end
    if (J(i) <= V1 - V1R) || (J(i) > V1 + V2R)
        new2 = (old1 + old2) / 2;
    else
        new2 = old2;
    end
    votes(I(i)) = new1;
    votes(J(i)) = new2;
end
error = (sum(votes(1 : V1) <= 0) + sum(votes(V1 + 1 : V1 + V2) >= 0)) / (V1 + V2);
duration = toc;
end