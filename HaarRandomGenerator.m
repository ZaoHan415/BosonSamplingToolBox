n = 100;
m = HaarUnitaryGenerator(n);
HaarDistributionValidation(m, n);


function out = HaarMatGen()
    out.generate = @HaarUnitaryGenerator;
    out.validate = @HaarDistributionValidation;
end

function retMat = HaarUnitaryGenerator(n)
    retMat = randn(n) + 1i*randn(n);
    [retMat, R] = qr(retMat);
    retMat = retMat * diag(sign(diag(R)));
end

function overlap = HaarDistributionValidation(haarMat, n)
    % entries r*exp(iq) of a Haar random matrix satisfy the distribution
    % 2(n-1)(1-r^2)^(n-2) r dr
    subplot(2,1,1);
    m = numel(haarMat);
    amplitudes = abs(reshape(haarMat, [1, m]));
    nbins = 100;
    h = histogram(amplitudes, [0:1/nbins:1], 'Normalization', 'probability');
    title('Distribution of amplitudes')
    hold on
    x = 0:1/nbins:1;
    haarDist = @(x) 2*(n-1)*(1-x.^2).^(n-2) .* x/nbins;
    plot(x, haarDist(x))
    legend('exper', 'theory')
    edges = h.BinEdges + 0.5/nbins;
    edges = edges(1:end-1);
    minVals = min(h.Values, haarDist(edges));
    overlap = sum(minVals);
    
    subplot(2,1,2);
    phases = angle(reshape(haarMat, [1, m]));
    phases = mod(phases, 2*pi);
    x = [0:2*pi/nbins:2*pi];
    histogram(phases, x, 'Normalization', 'probability');
    title('Distribution of phases')
    hold on
    fx1 = ones(1, length(x))/length(x);
    plot(x, fx1);
    legend('exper', 'theory')
end