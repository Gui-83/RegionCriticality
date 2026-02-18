function D = DiscreteBoundedPowerLawLikelihoodDerivative(n, l1, t, log_t, alpha)
    t_alpha = t.^(-alpha);
    l2 = sum(log_t.*t_alpha);
    l3 = sum(t_alpha);

    D = -l1 + n * l2 / l3;
end

