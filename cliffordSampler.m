% A = load('haarMat_10_1.mat');
% samples = cliffordSample(A.U(:, 1:2), 10)

function samples = cliffordSample(A, sampleSize)
    %A: n column of a m*m unitary matrix
    [m, n] = size(A);
    r = zeros(1, n);
    samples = zeros(n,sampleSize);
    for sampleCount=1:1:sampleSize
        if n>1
            A = A(:, randperm(n));
        end
        probs = abs(A(:, 1)).^2;
        r(1) = randsample(m, 1, true, probs);
        for k=2:1:n
            B = A(r(1:k-1), 1:k);
            pers = zeros(k, 1);
            for j=1:1:k
                colLabels = 1:k;
                colLabels(j) = [];
                pers(j) = permanent_Butler(B(:, colLabels));
            end
            probs = abs(A(:, 1:k) * pers).^2;
            r(k) = randsample(m, 1, true, probs);
        end
        samples(:, sampleCount) = sort(r)';
    end
end


% function w = sampleFromPMF(probs, n)
%     w = randsample(n, 10000, true, probs);
% end