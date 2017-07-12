function x = irwls(A, b, max_iter, epsilon)
% https://web.stanford.edu/~boyd/papers/pdf/rwl1.pdf
[~, p] = size(A);
w = ones(p, 1);


for i = 1:max_iter
  
  W = diag(w);
  
  cvx_begin quiet
    variable xhat(p,1)
    minimize(norm(W*xhat, 1))
    subject to
      A*xhat == b;
  cvx_end
  
  w = 1./(epsilon + abs(xhat));
  
end
x = xhat;
x(abs(x)<epsilon) = 0;


