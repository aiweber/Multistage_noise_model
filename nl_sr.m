function f = nl_sr(prs,x)
% function f = nl_sr(prs,x)
%
% Generates softplus function: 
%    f = prs(1)*log( 1 + exp(prs(2)*x + prs(3)) ) + prs(4)

f = zeros(size(x));

idx1 = prs(2)*x+prs(3)<500;
idx2 = prs(2)*x+prs(3)>=500;
f(idx1) = prs(1)*log(1+exp(prs(2)*x(idx1)+prs(3))) + prs(4);
f(idx2) = prs(1) * (prs(2)*x(idx2)+prs(3)) + prs(4);

