load("data");

%PLOT
plot(is_score, is_erro, 'ob', bp_score, bp_erro, 'or')
xlabel("Score");
ylabel("Classification error (%)");
legend("Algorithm 1", "Semi-supervised BP");