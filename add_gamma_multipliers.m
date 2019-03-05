function bem = add_gamma_multipliers(bem)


if length(bem.sigma) == 3
    sigma = [0 bem.sigma];
    bem.source_mult = 2./(sigma(2:end)+sigma(1:3));
    bem.field_mult = sigma(2:end)-sigma(1:3);
    num = (sigma(2:end)-sigma(1:3));
    den = (sigma(1:3)+sigma(2:end));
    bem.gamma = bsxfun(@rdivide,num,den');
elseif length(bem.sigma) == 1
    bem.source_mult = 2./bem.sigma;
    bem.field_mult = bem.sigma;
    bem.gamma = 1;
else
    error('meh')
end

