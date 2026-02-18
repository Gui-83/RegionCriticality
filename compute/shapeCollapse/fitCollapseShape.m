function gam = fitCollapseShape(T, shape, dicoStep)
    % compute the optimal gamma for the collapse shape by minimizing the variability of the common
    % shape acccros each difference lifetime
    % T = lifetime
    % shape = collapse shape (see docs/README.md for more information)
    % dicoStep = number of step to perform the dichotomy


    %n = size(shape, 2);
    %scalarProductMatrix = getScalarMatrix(shape);

    if isempty(shape)
        gam = NaN;
        return
    end

    left = 1;
    right = 4;
    for i = 1:dicoStep
        middle = (right+left)/2; % here middle is the candidate gamma
        %L_middle = Dloss(scalarProductMatrix, T, middle, n);
        f = scaleCollapseShape(shape, T, middle);
        L_middle = DvariabilityLoss(f, T);
        if L_middle < 0
            left = middle;
        else
            right = middle;
        end
    end
    gam = (left+right)/2;
end

function scalarProd = getScalarMatrix(shape)
    m = size(shape, 1);
    scalarProd = shape' * shape ./m;
end

function L = loss(scalarProductMatrix, T, gam, n)
    L = 0;
    for i = 1:n
        L1 = T(i)^(2-2*gam)*scalarProductMatrix(i, i);
        L2 = 1/n*T(i)^(1-gam)*sum(scalarProductMatrix(i,:).*T.^(1-gam));
        L3 = 1/n^2 * sum(scalarProductMatrix*transpose(T.^(1-gam))*T.^(1-gam), "all");

        norm_f_carre = T(i)^(2-2*gam)*scalarProductMatrix(i, i);
        g = L1 - 2*L2 + L3;
        L = L + g/norm_f_carre;
    end
end

function D = Dloss(scalarProductMatrix, T, gam, n)
    D = 0;
    for i = 1:n
        T_ratio = T(i)./T;
        L1 = 2/n*sum(scalarProductMatrix(i,:) .* log(T_ratio) .* T_ratio.^(gam-1));

        T_geometric_ratio = T(i) ./ sqrt(transpose(T)*T);
        L2 = 2/n^2*sum(scalarProductMatrix .* log(T_geometric_ratio) .* T_geometric_ratio.^(2*gam - 2), "all");
        norm_f_carre = T(i)^(1-gam)*scalarProductMatrix(i, i);
        dg = -L1 + L2;
        D = D + dg/norm_f_carre;
    end
end

function L = variabilityLoss(f)
    mu = mean(f, 2);
    variance = var(f, 0, 2);
    L = mean(variance ./ mu .^ 2);
end

function D = DvariabilityLoss(f, T)
    % this is the derivative of the variability at the point middle (in the main loop)
    m = size(f, 2);
    mu = mean(f, 2);
    f_ln_T = f.*log(T);
    mu_f_ln_T = mean(f_ln_T, 2);
    L1 = sum((f - mu).*(mu_f_ln_T - f_ln_T), 2);
    L2 = mu;
    L3 = mu_f_ln_T;
    L4 = var(f, 0, 2);

    L = (2/(m-1)*L1.*L2.^2 + 2*L2.*L3.*L4)./L2.^4;
    D = mean(L);

end