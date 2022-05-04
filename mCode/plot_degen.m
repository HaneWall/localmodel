clear all;
m = 9;
n = 1:1:((m-1)/2);

n_pu = 2*n;
n_pr = 1;
permutations = zeros(length(n),1);

for i=1:length(n)
    permutations(i) = multinomial_degen(m, n_pu(i), n_pr);
end


figure(1)
fig1_comps.fig = gcf;
p1 = stem(n, permutations/max(permutations));
set(p1, 'LineWidth', 1.5)
hold on
p2 = stem(n, permutations/(max(permutations)*m));
set(p2, 'LineWidth', 1.5)
set(gca, 'YScale', 'log')
xlim([0, m/2])
STANDARDIZE_FIGURE(fig1_comps)
