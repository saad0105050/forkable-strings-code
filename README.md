# forkable-strings-code
Compute the relative margin in polynomial time and space

# usage

# description
The programe takes the following as input:
* A sequence of positive integers n1, n2, ..., n_k, separated by space.
* A sequence of positive reals a1, a2, ..., a_m.

Let N = max{ n_1, ..., n_k } and let A = {alpha}
Let B(n, a) denote the binomial distribution with parameters n and a.
Let mu_x(y) be the relative margin

Our program essentially does the following:

for each n in { n_1, ..., n_k } do
	for each a in A do
		output Pr[mu_x(y) >= 0] where y ~ B(n, a)
	end for
end for
