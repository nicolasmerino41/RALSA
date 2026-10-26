ENV["OPPORTUNITY_STAGE"]="none"
include("opportunity_theory.jl")
using .OpportunityTheory, Random, Statistics, CSV, DataFrames
const T=OpportunityTheory
rng=MersenneTwister(917)
O=[[1,2,3],[2,3,4]];A=falses(2,4);A[1,1]=true;A[2,3]=true
s=T.SystemData("small_check",4,O,A,[3,3],[1,1],BitVector([true,true,false,false]))
p=(model="Site",a=-.5,h=1.)
states=collect(Iterators.product(ntuple(_->[-1.,1.],4)...))
qs=[T.logistic.(p.a .+ p.h.*collect(z)) for z in states]
ll=[sum(T.rowloglik(s,i,q[O[i]],s.keep) for i in 1:2) for q in qs]
w=exp.(ll.-T.lse(ll));exact=zeros(2,4)
for (q,weight) in zip(qs,w),i in 1:2
    test=filter(j->!s.keep[j],O[i]);positive=!any(s.A[i,s.keep])
    norm=positive ? -expm1(sum(log1p.(-q[test]))) : 1.
    for j in test;exact[i,j]+=weight*q[j]/norm;end
end
draws,chains=T.withinpredict(rng,s,p)
estimated=reduce(+,Float64.(x) for x in draws)/length(draws)
error=maximum(abs.(estimated[:,.!s.keep].-exact[:,.!s.keep]))
@assert error<.1
for model in T.MODELS
    samples=[T.simulate(rng,s,model,-.5,1.) for _ in 1:10000]
    @assert all(all(vec(sum(x;dims=2)).>=1) for x in samples)
    empirical=mean(sum(x) for x in samples)
    expected=T.meanactivity(s,model,-.5,1.,zeros(4,1))
    @assert abs(empirical-expected)<.06
end
CSV.write(joinpath(T.DEST,"conditional_sampler_verification.csv"),DataFrame(check=["Site posterior against exhaustive 16-state integration","Positive-row sampler mean against exact expectation"],max_error=[error,NaN],passed=[true,true]))
println("Conditional sampler verified; largest probability discrepancy = ",error)
