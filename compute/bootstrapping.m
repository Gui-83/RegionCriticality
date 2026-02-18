function [BS, BT, BST] = bootstrapping(S, T, ST, samples, k)
    %bootstrap S, T and ST k times withe 'samples' samples each
    BS = zeros(samples, k);
    BT = zeros(samples, k);
    BST = {};
    n = size(S, 1);
    for i = 1:k
        indices = randi([1 n], 1, samples);
        BS(:,i) = S(indices);
        BT(:,i) = T(indices);
        BST{i} = ST(indices);
    end
end

