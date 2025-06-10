%This function generates a matrix A which is upper triangular and stored in
%sparse form. It corresponds to the ajacency matrix of a stochastic block
%the sizes of the communities are stored in the column vector n and the
%edge probabilities are stored in the matrix P.

function [A, I, J] = sbm(n, P)
    %Matrix entries with ones
    I = zeros(1, round(1.1 * n' * P * n / 2));
    J = zeros(1, round(1.1 * n' * P * n / 2));
    %Starting indexes of communities
    s = cumsum([1; n(1 : end)]);
    %Generating matrix
    k = 0;
    for a = 1 : length(n)
        for b = a : length(n)
            for i = s(a) : s(a + 1) - 1
                j = max([i + 1 s(b)]) + geornd(P(a, b));
                while j < s(b + 1)
                    k = k + 1;
                    I(k) = i;
                    J(k) = j;
                    j = j + 1 + geornd(P(a, b));
                end
            end
        end
    end
    %Creating sparse matrix
    I = I(1 : k);
    J = J(1 : k);
    A = sparse(I, J, ones(1, k), sum(n), sum(n));
end