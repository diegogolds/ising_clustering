load("log_3");
e = 100 * e;

%AVERAGES AND VARIANCES
dav = zeros(10, 2);
eav = zeros(10, 2);
iav = zeros(10, 2);
dsd = zeros(10, 2);
esd = zeros(10, 2);
isd = zeros(10, 2);
for k = 1 : 10
    for l = 1 : 8
        dav(k, l) = sum(d(k, l, :)) / rounds;
        eav(k, l) = sum(e(k, l, :)) / rounds;
        iav(k, l) = sum(i(k, l, :)) / rounds;
        dsd(k, l) = sqrt(sum((d(k, l, :) - dav(k, l)) .^ 2) / rounds);
        esd(k, l) = sqrt(sum((e(k, l, :) - eav(k, l)) .^ 2) / rounds);
        isd(k, l) = sqrt(sum((i(k, l, :) - iav(k, l)) .^ 2) / rounds);
    end
end

%PLOTS
figure;
plot(eta, dav(:, 1), 'o', eta, dav(:, 2), '^', eta, dav(:, 3), 'square', eta, dav(:, 4), 'diamond', eta, dav(:, 5), 'pentagram', eta, dav(:, 6), 'x', eta, dav(:, 7), '+', eta, dav(:, 8), '*');
title("Logarithmic degrees with n = 5000, a = 3 and b = 1");
legend('\alpha = 10 and \beta = \infty', 'Asynchronous Consensus', 'Synchronous Consensus', 'Gossiping', 'Standard Laplacian', 'Normalized Laplacian', 'Page Rank', 'Poisson Learning');
xlabel('\eta');
ylabel('Running time (sec)');
figure;
plot(eta, eav(:, 1), 'o', eta, eav(:, 2), '^', eta, eav(:, 3), 'square', eta, eav(:, 4), 'diamond', eta, eav(:, 5), 'pentagram', eta, eav(:, 6), 'x', eta, eav(:, 7), '+', eta, eav(:, 8), '*');
title("Logarithmic degrees with n = 5000, a = 3 and b = 1");
legend('\alpha = 10 and \beta = \infty', 'Asynchronous Consensus', 'Synchronous Consensus', 'Gossiping', 'Standard Laplacian', 'Normalized Laplacian', 'Page Rank', 'Poisson Learning');
xlabel('\eta');
ylabel('Relative error (%)');