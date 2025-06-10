clear all;
load("n5000a5eta01loglog")

%Experiment selection
X = Xzri;
Y = Yzri;
%Discarding bad realizations
X = [X zeros(rounds, 1)];
Y = [Y zeros(rounds, 1)];
X(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
Y(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
samples = sum(X(:, end));
disp("Fail rate = " + 100 * (rounds - samples) / rounds + "%");
X = X(X(:, end) == 1, 1 : end - 1) / n;
Y = Y(Y(:, end) == 1, 1 : end - 1) / n;
%Sample means
Xav0 = sum(X) / samples;
Yav0 = sum(Y) / samples;
%Sample variances
Xvr = sum((X - Xav0) .^ 2) / samples;
Yvr = sum((Y - Yav0) .^ 2) / samples;
%Confidence level
c = norminv(0.975);
%Confidence intervals
Xup0 = Xav0 + c * sqrt(Xvr / samples);
Xlw0 = Xav0 - c * sqrt(Xvr / samples);
Yup0 = Yav0 + c * sqrt(Yvr / samples);
Ylw0 = Yav0 - c * sqrt(Yvr / samples);

load("n5000a5eta01sqrtlog")

%Experiment selection
X = Xzri;
Y = Yzri;
%Discarding bad realizations
X = [X zeros(rounds, 1)];
Y = [Y zeros(rounds, 1)];
X(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
Y(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
samples = sum(X(:, end));
disp("Fail rate = " + 100 * (rounds - samples) / rounds + "%");
X = X(X(:, end) == 1, 1 : end - 1) / n;
Y = Y(Y(:, end) == 1, 1 : end - 1) / n;
%Sample means
Xav1 = sum(X) / samples;
Yav1 = sum(Y) / samples;
%Sample variances
Xvr = sum((X - Xav1) .^ 2) / samples;
Yvr = sum((Y - Yav1) .^ 2) / samples;
%Confidence level
c = norminv(0.975);
%Confidence intervals
Xup1 = Xav1 + c * sqrt(Xvr / samples);
Xlw1 = Xav1 - c * sqrt(Xvr / samples);
Yup1 = Yav1 + c * sqrt(Yvr / samples);
Ylw1 = Yav1 - c * sqrt(Yvr / samples);

load("n5000a5eta01log")

%Experiment selection
X = Xzri;
Y = Yzri;
%Discarding bad realizations
X = [X zeros(rounds, 1)];
Y = [Y zeros(rounds, 1)];
X(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
Y(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
samples = sum(X(:, end));
disp("Fail rate = " + 100 * (rounds - samples) / rounds + "%");
X = X(X(:, end) == 1, 1 : end - 1) / n;
Y = Y(Y(:, end) == 1, 1 : end - 1) / n;
%Sample means
Xav2 = sum(X) / samples;
Yav2 = sum(Y) / samples;
%Sample variances
Xvr = sum((X - Xav2) .^ 2) / samples;
Yvr = sum((Y - Yav2) .^ 2) / samples;
%Confidence level
c = norminv(0.975);
%Confidence intervals
Xup2 = Xav2 + c * sqrt(Xvr / samples);
Xlw2 = Xav2 - c * sqrt(Xvr / samples);
Yup2 = Yav2 + c * sqrt(Yvr / samples);
Ylw2 = Yav2 - c * sqrt(Yvr / samples);

load("n5000a5eta01sqrt")

%Experiment selection
X = Xzri;
Y = Yzri;
%Discarding bad realizations
X = [X zeros(rounds, 1)];
Y = [Y zeros(rounds, 1)];
X(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
Y(:, end) = X(:, end - 1) > 0.1 & Y(:, end - 1) < -0.1;
samples = sum(X(:, end));
disp("Fail rate = " + 100 * (rounds - samples) / rounds + "%");
X = X(X(:, end) == 1, 1 : end - 1) / n;
Y = Y(Y(:, end) == 1, 1 : end - 1) / n;
%Sample means
Xav3 = sum(X) / samples;
Yav3 = sum(Y) / samples;
%Sample variances
Xvr = sum((X - Xav3) .^ 2) / samples;
Yvr = sum((Y - Yav3) .^ 2) / samples;
%Confidence level
c = norminv(0.975);
%Confidence intervals
Xup3 = Xav3 + c * sqrt(Xvr / samples);
Xlw3 = Xav3 - c * sqrt(Xvr / samples);
Yup3 = Yav3 + c * sqrt(Yvr / samples);
Ylw3 = Yav3 - c * sqrt(Yvr / samples);

%Plot
figure;
plot(t, Xav0, 'y', t, Xav1, 'b', t, Xav2, 'r', t, Xav3, 'g', t, 1 + (eta - 1)* exp(-t), ':k', t, -1 - (eta - 1)* exp(-t), ':k', t, Xup0, ':y', t, Xlw0, ':y', t, Yav0, 'y', t, Yup0, ':y', t, Ylw0, ':y', t, Xup1, ':b', t, Xlw1, ':b', t, Yav1, 'b', t, Yup1, ':b', t, Ylw1, ':b', t, Xup2, ':r', t, Xlw2, ':r', t, Yav2, 'r', t, Yup2, ':r', t, Ylw2, ':r', t, Xup3, ':g', t, Xlw3, ':g', t, Yav3, 'g', t, Yup3, ':g', t, Ylw3, ':g');
title("Magnetization over time for \alpha = 0 and \beta = \infty");
legend("\lambda_n = loglog(n)", "\lambda_n = [log(n)]^{1/2}", "\lambda_n = log(n)", "\lambda_n = n^{1/2}", "Mean-field limit");
ylabel('Magnetization');
xlabel('t');

