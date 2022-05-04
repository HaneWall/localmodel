function [degen] = multinomial(m, n_pr_d, n_pr_u, n_pu_d, n_pu_u)
degen = factorial(m)/(factorial(n_pr_u) * factorial(n_pr_d) * factorial(n_pu_d) * factorial(n_pu_u));
end