# Exact probability that all k supporting sites are among m removed sites.
loss(k, N, m) = k > m ? 0.0 : prod(((m-j)/(N-j) for j in 0:k-1); init=1.0)
removed_count(N, f) = N - max(1, round(Int, N * (1-f)))
