function gam = fitShapeCollapse(profile,dicoStep)
% compute the optimal gamma from shape collapse by minimizing the variability of the common shape acccros each difference lifetime
% profile : cell array, each element is an avalanche profile (column vector)
% dicoStep : number of step to perform the dichotomy

arguments
  profile
  dicoStep (1,1) = 15
end

if isempty(profile)
  gam = NaN;
  return
end

% average profile for each unique lifetime
size_t_avrg = AvalAverageSizeTimeDependent(profile);
% interpolate every shape to a common rescaled time
[~,shape,T] = transformCollapseShape(size_t_avrg);
T = reshape(T,1,[]);

% dichotomy
left = 1;
right = 4;
for i = 1 : dicoStep
  middle = (right+left) / 2; % candidate gamma
  f = scaleCollapseShape(shape, T, middle);
  L_middle = DvariabilityLoss(f,T);
  if L_middle < 0
    left = middle;
  else
    right = middle;
  end
end
gam = (left+right) / 2;

end

% --- helper functions ---

function scale_shape = scaleCollapseShape(shape, T, gam)
% scale the shape according to eq. 49 of Yang Tian Theoretical foundations of studying ...
  scale_shape = shape.*T.^(1-gam);
end

function D = DvariabilityLoss(f, T)
    % derivative of the variability at the middle point (in the main loop)
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