load("10_6_nodes_eta_06")

error_percen = 100 * target_error;

%Averages
is_mean_iter = mean(is_iter, 2);
is_mean_flip = mean(is_flip, 2);
bp_mean_iter = mean(bp_iter, 2);

%Scores
is_score = is_mean_iter - is_mean_flip + (3 + (a + b) * lambda(n)) * is_mean_flip;
bp_score = (2 * (1 - eta)^2 * (a + b) * n * lambda(n) + 2 * (1 - eta) * n) * bp_mean_iter;

%Comparison
figure;
plot(error_percen, is_score, 'ob', error_percen, bp_score, 'or');
legend("Ising", "Belief Propagation");
xlabel("Target error (%)");
ylabel("Score");

disp(bp_score ./ is_score);
