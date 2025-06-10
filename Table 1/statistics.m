load("log_10");
e = 100 * e;

%AVERAGES AND VARIANCES
dav = zeros(10, 2);
eav = zeros(10, 2);
iav = zeros(10, 2);
dsd = zeros(10, 2);
esd = zeros(10, 2);
isd = zeros(10, 2);
for k = 1 : 10
    for l = 1 : 4
        dav(k, l) = sum(d(k, l, :)) / rounds;
        eav(k, l) = sum(e(k, l, :)) / rounds;
        iav(k, l) = sum(i(k, l, :)) / rounds;
        dsd(k, l) = sqrt(sum((d(k, l, :) - dav(k, l)) .^ 2) / rounds);
        esd(k, l) = sqrt(sum((e(k, l, :) - eav(k, l)) .^ 2) / rounds);
        isd(k, l) = sqrt(sum((i(k, l, :) - iav(k, l)) .^ 2) / rounds);
    end
end