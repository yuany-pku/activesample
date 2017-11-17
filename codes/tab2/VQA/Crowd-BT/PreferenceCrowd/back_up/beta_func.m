function y = beta_func(z,w)
if nargin<2,
  error(message('MATLAB:beta:NotEnoughInputs'));
end
y = exp(gammaln(z)+gammaln(w)-gammaln(z+w));
end
