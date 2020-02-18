function f = nl_sr_reparam(prs,x)
% function f = nl_sr_reparam(prs,x)
%
% Generates softplus function, adjusted to match reparameterization
% in optimization code:
%    f = prs(1)*log( 1 + exp(prs(2)/prs(1)*x(idx1) + prs(3)/prs(1)) ) + prs(4)

if length(prs)<4
    prs(4) = 0;
end

f = zeros(size(x));

idx1 = prs(2)/prs(1)*x+prs(3)/prs(1)<500;
idx2 = prs(2)/prs(1)*x+prs(3)/prs(1)>=500;

f(idx1) = prs(1)*log(1+exp(prs(2)/prs(1)*x(idx1)+prs(3)/prs(1))) + prs(4);
f(idx2) = prs(1) * (prs(2)/prs(1)*x(idx2)+prs(3)/prs(1)) + prs(4);

