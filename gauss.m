function y = gauss(x, mu, sigma)
% y = gauss(x, mu, sigma)

y = exp(-(x-mu).^2/(2*sigma^2));
