function logcombinations(C)
 [vcat(0.,cumsum([log((c-k+1)/k) for k in 1:c])) for c in C]
end
function sweep!(rng,K,C,active,lc)
 P=length(K)
 for _ in 1:P
  a=rand(rng,eachindex(active));i=active[a];j=rand(rng,1:P);i==j && continue
  t=K[i]+K[j]
  if K[j]==0
   if t<=C[j]
    z=lc[j][t+1]-lc[i][t+1]
    if rand(rng)<1/(1+exp(-z));K[i]=0;K[j]=t;active[a]=j;end
   end
  else
   lo=max(1,t-C[j]);hi=min(C[i],t-1)
   lo==hi && continue
   mx=-Inf
   for k in lo:hi;mx=max(mx,lc[i][k+1]+lc[j][t-k+1]);end
   total=0.
   for k in lo:hi;total+=exp(lc[i][k+1]+lc[j][t-k+1]-mx);end
   u=rand(rng)*total;v=0.;chosen=hi
   for k in lo:hi
    v+=exp(lc[i][k+1]+lc[j][t-k+1]-mx)
    if v>=u;chosen=k;break;end
   end
   K[i]=chosen;K[j]=t-chosen
  end
 end
end
function diagnostics(x)
 n,c=size(x);half=n÷2;y=hcat(x[1:half,:],x[end-half+1:end,:]);n=size(y,1)
 W=mean(var(y;dims=1));between=n*var(vec(mean(y;dims=1)))
 rhat=W>0 ? sqrt(((n-1)*W/n+between/n)/W) : 1.
 ess=0.
 for j in 1:c
  v=x[:,j].-mean(x[:,j]);den=sum(abs2,v);tau=1.
  if den>0
   for lag in 1:2:min(length(v)÷2,500)-1
    rho=sum(v[1:end-lag].*v[1+lag:end])/den+sum(v[1:end-lag-1].*v[2+lag:end])/den
    rho<=0 && break;tau+=2rho
   end
  end
  ess+=length(v)/tau
 end
 (rhat=rhat,ess=ess,mcse=std(vec(x))/sqrt(ess))
end
function sampler_test()
 C=[2,3,4,2];M=6;L=3;lc=logcombinations(C)
 states=[collect(k) for k in Iterators.product((0:c for c in C)...) if sum(k)==M && count(>(0),k)==L]
 w=[prod(binomial(c,k) for (c,k) in zip(C,K)) for K in states];p=w/sum(w)
 lookup=Dict(Tuple(k)=>i for (i,k) in enumerate(states));Q=zeros(length(states),length(states))
 for (a,K) in enumerate(states),i in findall(>(0),K),j in eachindex(C)
  rate=1/(L*length(C))
  if i==j;Q[a,a]+=rate;continue;end
  total=K[i]+K[j];num=(K[i]>0)+(K[j]>0)
  options=[k for k in max(0,total-C[j]):min(C[i],total) if (k>0)+(total-k>0)==num]
  weights=[binomial(C[i],k)*binomial(C[j],total-k) for k in options];weights=weights/sum(weights)
  for (k,prob) in zip(options,weights)
   v=copy(K);v[i]=k;v[j]=total-k;Q[a,lookup[Tuple(v)]]+=rate*prob
  end
 end
 err=maximum(abs.(p.*Q-transpose(p.*Q)));@assert err<1e-12 && maximum(abs.(sum(Q;dims=2).-1))<1e-12
 CSV.write(joinpath(OUT,"sampler_check.csv"),DataFrame(states=[length(states)],detailed_balance_error=[err],passed=[true]))
end
