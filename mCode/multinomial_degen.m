function [degeneracy] = multinomial_degen(m, n_pu, n_pr)
    %m needs has to be odd
    degeneracy = factorial(m)/(factorial(n_pr)*factorial((m + n_pu - n_pr)/2)*factorial((m - n_pu - n_pr)/2));
end
