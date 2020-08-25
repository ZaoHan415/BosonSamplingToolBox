% using a naive two-step-procedure to assign each output configuration a
% unique integer
% e.g. consider n=3,m=4 case
n = 2;
m = 10;
sampleNum = 100000;
A = load('haarMat_10_1.mat');
U = A.U(:, 1:2);


samples = cliffordSampler(U, sampleNum);
total = nchoosek(n+m-1, n);
binCounts = zeros(1, total);
[probs, orderList] = getOrder(n, m, U);
for k=1:1:sampleNum
    t = orderList(bitSum(samples(:, k), m));
    binCounts(t) = binCounts(t) + 1;
end
subplot(2,1,1);
histogram('BinEdges',0:total,'BinCounts',binCounts)
title("distribution from Clifford-Clifford sampling")
subplot(2,1,2);
histogram('BinEdges',0:total,'BinCounts',probs)
title("distribution from direct calculation")

function [probs, orderList] = getOrder(n, m, U)
    total = nchoosek(n+m-1, n);
    probs = zeros(1, total);
    seq = zeros(1, total);
    config = ones(1, n);
    for k=1:1:total
        seq(k) = bitSum(config, m);
        probs(k) = abs(permanent_Butler(U(config, :))).^2 / multipli(config, m, n);
        if k == total
            break
        end
        config = update(config, m);
    end
    orderList = containers.Map(seq, 1:1:total);
end

function result = update(config, m)
    if length(config) == 1
        result = config + 1;
    end
    t = config(1) + 1;
    result = config;
    if t > m
        result(2:end) = update(config(2:end), m);
        result(1) = result(2);
    else
        result(1) = result(1) + 1;
    end
end

function order = bitSum(lst, m)
    order = sum(m.^(lst-1));
end

function t = multipli(config, m, n)
    lst = zeros(1, m);
    for i=1:1:n
        lst(config(i)) = lst(config(i)) + 1;
    end
    t = prod(factorial(lst));
end
    