# Local-support distribution analysis (Julia 1.12; CSV, DataFrames, CairoMakie, SHA, Statistics, Random)
# Location: <repo>/scripts/local_support_distributions.jl
# Usage:
#   julia --project=. scripts/local_support_distributions.jl
#   julia --project=. scripts/local_support_distributions.jl --scenario inclusive --replicates 200
#   julia --project=. scripts/local_support_distributions.jl --subsamples 0 --replicates 0
#   julia --project=. scripts/local_support_distributions.jl --root /path/to/repo
# Primary models use ALL co-occurring pairs (K >= 0) and their fixed observed opportunities n.
# The separate descriptive distribution fits use only K >= 1, and are truncated to 1:N.

module LocalSupportDistributions

using CSV, DataFrames, Statistics, Random, SHA, CairoMakie

const DATA = ["Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
              "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"]
const COLORS = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100"]
const DEPS_DESCRIPTION = "Only CSV, DataFrames, CairoMakie, Random, Statistics and SHA are used."

struct Settings
    root::String
    scenario::String
    replicates::Int
    subsamples::Int
    seed::Int
end

function parse_args(args)
    root = normpath(joinpath(@__DIR__, ".."))
    scenario = "confirmed"
    replicates = 100
    subsamples = 40
    seed = 20261014
    i = 1
    while i <= length(args)
        key = args[i]
        key in ("--root", "--scenario", "--replicates", "--subsamples", "--seed") ||
            error("Unknown argument: $key")
        i < length(args) || error("Missing value for $key")
        value = args[i+1]
        if key == "--root"
            root = abspath(value)
        elseif key == "--scenario"
            scenario = value
        elseif key == "--replicates"
            replicates = parse(Int, value)
        elseif key == "--subsamples"
            subsamples = parse(Int, value)
        else
            seed = parse(Int, value)
        end
        i += 2
    end
    scenario in ("confirmed", "inclusive") || error("--scenario must be confirmed or inclusive")
    replicates >= 0 && subsamples >= 0 || error("Simulation counts cannot be negative")
    return Settings(root, scenario, replicates, subsamples, seed)
end

file_sha(path) = bytes2hex(sha256(read(path)))
function store(path, table)
    mkpath(dirname(path))
    CSV.write(path, table)
end

# Pure-Julia positive-argument Lanczos log-gamma (avoids new Project.toml dependencies).
# Valid for x>0; tested numerically against known factorial identities below.
const LCOEF = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
               771.32342877765313, -176.61502916214060, 12.507343278686905,
               -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
function loggamma_pos(x::Real)
    z = Float64(x)
    z > 0 || throw(DomainError(x, "loggamma_pos expects positive arguments"))
    if z < 0.5
        return log(pi) - log(sinpi(z)) - loggamma_pos(1.0-z)
    end
    zz = z - 1.0
    a = LCOEF[1]
    @inbounds for j in 2:length(LCOEF)
        a += LCOEF[j] / (zz + (j-1))
    end
    t = zz + 7.5
    return 0.5*log(2pi) + (zz+0.5)*log(t) - t + log(a)
end
logchoose(n::Int, k::Int) = loggamma_pos(n+1) - loggamma_pos(k+1) - loggamma_pos(n-k+1)
logbeta(a,b) = loggamma_pos(a) + loggamma_pos(b) - loggamma_pos(a+b)
sigmoid(x) = x >= 0 ? 1/(1+exp(-x)) : exp(x)/(1+exp(x))
logit(p) = log(p/(1-p))
clamp_prob(p) = clamp(Float64(p), 1e-9, 1-1e-9)

struct Cell
    n::Int
    k::Int
    count::Int
    logcomb::Float64
end

function cells_from(v::DataFrame)
    cc = combine(groupby(v, [:n, :K]), nrow => :count)
    return [Cell(Int(r.n), Int(r.K), Int(r.count), logchoose(Int(r.n), Int(r.K))) for r in eachrow(cc)]
end

function merge_cells(vectors)
    counts=Dict{Tuple{Int,Int},Int}()
    for cells in vectors, c in cells
        key=(c.n,c.k)
        counts[key]=get(counts,key,0)+c.count
    end
    return [Cell(n,k,count,logchoose(n,k)) for ((n,k),count) in counts]
end

function load_or_build_pairs(input, out)
    pairspath=joinpath(input,"pairs.csv")
    if isfile(pairspath)
        return CSV.read(pairspath,DataFrame), true
    end
    @warn "pairs.csv not present for scenario: reconstructing from local co-occurrence and interaction tables" input
    rows=NamedTuple[]
    for ds in DATA
        base=joinpath(input,"inputs")
        co=CSV.read(joinpath(base,ds*"_cooccurrences.csv"),DataFrame;types=String)
        ints=CSV.read(joinpath(base,ds*"_interactions.csv"),DataFrame;types=String)
        occ=CSV.read(joinpath(base,ds*"_occupancy.csv"),DataFrame;types=String)
        sites=length(unique(occ.site))
        n=Dict{Tuple{String,String},Int}()
        k=Dict{Tuple{String,String},Int}()
        for r in eachrow(co)
            key=(r.consumer,r.resource)
            n[key]=get(n,key,0)+1
        end
        for r in eachrow(ints)
            key=(r.consumer,r.resource)
            k[key]=get(k,key,0)+1
        end
        for ((consumer,resource),opportunities) in n
            push!(rows,(dataset=ds,consumer=consumer,resource=resource,
                        n=opportunities,K=get(k,(consumer,resource),0),N=sites))
        end
    end
    result=DataFrame(rows)
    sort!(result,[:dataset,:consumer,:resource])
    store(joinpath(out,"pairs_rebuilt_from_inputs.csv"),result)
    return result,false
end

function verify_data(input, pairs)
    required = [:dataset, :consumer, :resource, :K, :n, :N]
    all(x -> x in propertynames(pairs), required) || error("pairs.csv is missing one of $required")
    any(ismissing, Matrix(pairs[:, required])) && error("pairs.csv has missing values")
    @assert Set(String.(unique(pairs.dataset))) == Set(DATA) "Unexpected dataset list in pairs.csv"
    @assert isapprox(loggamma_pos(1), 0.0; atol=1e-12)
    @assert isapprox(loggamma_pos(6), log(120.0); atol=1e-11)
    for ds in DATA
        v = pairs[pairs.dataset .== ds, :]
        nrow(v) > 0 || error("Dataset $ds has no pairs")
        length(unique(v.N)) == 1 || error("Inconsistent site totals for $ds")
        nrow(unique(v[:, [:consumer, :resource]])) == nrow(v) || error("Duplicate pairs in $ds")
        all((v.K .>= 0) .& (v.n .>= 1) .& (v.K .<= v.n) .& (v.n .<= v.N)) ||
            error("Invalid 0<=K<=n<=N for $ds")
    end
    # Verify source hashes, as done in EmpiricalSupportScatter.
    versionpath=joinpath(input,"input_versions.csv")
    if isfile(versionpath)
        ver=CSV.read(versionpath,DataFrame)
        for r in eachrow(ver)
            file=joinpath(input,"inputs",string(r.dataset,"_",r.kind,".csv"))
            isfile(file) || error("Missing input: $file")
            file_sha(file)==string(r.sha256) || error("SHA-256 mismatch: $file")
        end
    else
        @warn "No input_versions.csv in this scenario; computing provenance hashes without comparison to a frozen manifest" input
    end
    # Reconstruct n and K directly from site-level presence/interaction records.
    # Every input row is a binary per-site event; duplicates would be an error.
    for ds in DATA
        co = CSV.read(joinpath(input, "inputs", ds*"_cooccurrences.csv"), DataFrame; types=String)
        ints = CSV.read(joinpath(input, "inputs", ds*"_interactions.csv"), DataFrame; types=String)
        occ = CSV.read(joinpath(input, "inputs", ds*"_occupancy.csv"), DataFrame; types=String)
        for tbl in (co, ints)
            all(x -> x in propertynames(tbl), [:site,:consumer,:resource]) || error("Bad site-level schema: $ds")
            nrow(unique(tbl)) == nrow(tbl) || error("Duplicate site pairs: $ds")
        end
        N = length(unique(occ.site))
        v = pairs[pairs.dataset .== ds, :]
        N == Int(v.N[1]) || error("N differs from occupancy for $ds: $N != $(v.N[1])")
        nc = Dict{Tuple{String,String},Int}()
        kc = Dict{Tuple{String,String},Int}()
        for r in eachrow(co)
            key = (r.consumer,r.resource)
            nc[key] = get(nc,key,0)+1
        end
        for r in eachrow(ints)
            key = (r.consumer,r.resource)
            kc[key] = get(kc,key,0)+1
        end
        length(nc) == nrow(v) || error("Pair count disagrees for $ds")
        for r in eachrow(v)
            key = (string(r.consumer),string(r.resource))
            Int(r.n) == get(nc,key,-1) || error("Opportunity mismatch for $ds $key")
            Int(r.K) == get(kc,key,0) || error("Support mismatch for $ds $key")
        end
        println("  verified ", ds, ": ", nrow(v), " co-occurring pairs, ", count(v.K .> 0), " positive links")
    end
end

function basic_distributions(data::Dict{String,DataFrame}, out)
    pmf = NamedTuple[]; opportunities = NamedTuple[]; cond = NamedTuple[]
    conditional_fraction = NamedTuple[]; desc = NamedTuple[]
    binlabel(n) = n==1 ? "1" : n<=3 ? "2-3" : n<=7 ? "4-7" : n<=15 ? "8-15" : n<=31 ? "16-31" : n<=63 ? "32-63" : "64+"
    bins = ("1", "2-3", "4-7", "8-15", "16-31", "32-63", "64+")
    for ds in DATA
        v=data[ds]; N=Int(v.N[1]); allcount=nrow(v); positive=count(v.K .> 0)
        positive > 0 || error("No positive links in $ds: cannot fit positive-support distributions")
        cs=zeros(Int,N+1); ns=zeros(Int,N+1)
        for r in eachrow(v)
            cs[Int(r.K)+1]+=1; ns[Int(r.n)+1]+=1
        end
        for k in 0:N
            for onlypositive in (false,true)
                if !onlypositive || k>=1
                    denominator = onlypositive ? positive : allcount
                    freq = onlypositive ? (k==0 ? 0 : cs[k+1]) : cs[k+1]
                    tail = sum(cs[max(k, onlypositive ? 1 : 0)+1:end]) / denominator
                    push!(pmf, (dataset=ds,positive_only=onlypositive,K=k,
                                count=freq,probability=freq/denominator,ccdf=tail,
                                denominator=denominator))
                end
            end
            if k>=1
                push!(opportunities,(dataset=ds,n=k,count=ns[k+1],probability=ns[k+1]/allcount,
                                     ccdf=sum(ns[k+1:end])/allcount))
            end
        end
        for b in bins
            inds = [i for i in 1:nrow(v) if binlabel(Int(v.n[i]))==b]
            isempty(inds) && continue
            for posonly in (false,true)
                selected=[i for i in inds if !posonly || Int(v.K[i])>0]
                isempty(selected) && continue
                cc=Dict{Int,Int}()
                for i in selected
                    k=Int(v.K[i]); cc[k]=get(cc,k,0)+1
                end
                for k in sort(collect(keys(cc)))
                    push!(cond,(dataset=ds,n_bin=b,positive_only=posonly,K=k,
                                count=cc[k],probability=cc[k]/length(selected),
                                n_pairs=length(selected)))
                end
                fraction_counts=zeros(Int,5)
                for i in selected
                    f=Int(v.K[i])/Int(v.n[i])
                    j=f==0 ? 1 : f<=0.25 ? 2 : f<=0.50 ? 3 : f<=0.75 ? 4 : 5
                    fraction_counts[j]+=1
                end
                labels=("0", "(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")
                for j in 1:5
                    push!(conditional_fraction,(dataset=ds,n_bin=b,positive_only=posonly,
                          fraction_bin=labels[j],probability=fraction_counts[j]/length(selected),
                          count=fraction_counts[j],n_pairs=length(selected)))
                end
            end
        end
        pos = v[v.K .> 0,:]
        push!(desc,(dataset=ds,sites=N,cooccurring_pairs=allcount,regional_links=positive,
                    zero_fraction=cs[1]/allcount,singleton_fraction_positive=cs[2]/positive,
                    doubleton_or_less_positive=sum(cs[2:min(N+1,3)])/positive,
                    mean_K_positive=mean(pos.K),median_K_positive=median(pos.K),
                    mean_n_all=mean(v.n),mean_fraction_all=mean(v.K ./ v.n),
                    mean_fraction_positive=mean(pos.K ./ pos.n),
                    total_support=sum(v.K),total_opportunities=sum(v.n)))
    end
    # Pairwise similarities of marginal K distributions (deliberately NOT opportunity-adjusted).
    pairwise=NamedTuple[]
    full=DataFrame(pmf)
    for a in DATA, b in DATA
        ta=full[(full.dataset.==a).&(.!full.positive_only),:]
        tb=full[(full.dataset.==b).&(.!full.positive_only),:]
        maxN=max(maximum(ta.K),maximum(tb.K))
        pa=zeros(Float64,maxN+1); pb=zeros(Float64,maxN+1)
        pa[Int.(ta.K).+1]=ta.probability; pb[Int.(tb.K).+1]=tb.probability
        push!(pairwise,(dataset_A=a,dataset_B=b,JS_marginal_K=js_divergence(pa,pb)))
    end
    store(joinpath(out,"01_empirical_distributions","pairwise_marginal_JS.csv"),DataFrame(pairwise))
    store(joinpath(out,"01_empirical_distributions","pmf_ccdf.csv"),DataFrame(pmf))
    store(joinpath(out,"01_empirical_distributions","opportunity_pmf_ccdf.csv"),DataFrame(opportunities))
    store(joinpath(out,"01_empirical_distributions","conditional_K_by_n_bin.csv"),DataFrame(cond))
    store(joinpath(out,"01_empirical_distributions","conditional_fraction_by_n_bin.csv"),DataFrame(conditional_fraction))
    store(joinpath(out,"01_empirical_distributions","dataset_summary.csv"),DataFrame(desc))
    return DataFrame(pmf),DataFrame(opportunities),DataFrame(desc)
end

# Bound-constrained 1D golden section maximization.
function golden_max(f, low, high; iterations=100)
    g=(sqrt(5.0)-1)/2
    a=Float64(low); b=Float64(high)
    c=b-g*(b-a); d=a+g*(b-a)
    fc=f(c); fd=f(d)
    for _ in 1:iterations
        if fc>fd
            b=d; d=c; fd=fc; c=b-g*(b-a); fc=f(c)
        else
            a=c; c=d; fc=fd; d=a+g*(b-a); fd=f(d)
        end
    end
    x=(a+b)/2
    return (x,f(x))
end

# Small, dependency-free, two-dimensional Nelder-Mead. Parameters optimized on unconstrained scales.
function nm_min(f, x0; iterations=280, tol=1e-8)
    pts=[copy(Float64.(x0)),Float64.(x0).+[0.35,0.0],Float64.(x0).+[0.0,0.35]]
    vals=[f(x) for x in pts]
    for _ in 1:iterations
        order=sortperm(vals); pts=pts[order]; vals=vals[order]
        if maximum(norm2(pts[j].-pts[1]) for j in 2:3)<tol && abs(vals[3]-vals[1])<tol
            break
        end
        centroid=(pts[1].+pts[2])./2
        reflect=centroid.+(centroid.-pts[3]); fr=f(reflect)
        if fr<vals[1]
            expand=centroid.+2 .* (reflect.-centroid); fe=f(expand)
            pts[3],vals[3]=fe<fr ? (expand,fe) : (reflect,fr)
        elseif fr<vals[2]
            pts[3]=reflect; vals[3]=fr
        else
            contract=fr<vals[3] ? centroid.+0.5.*(reflect.-centroid) : centroid.+0.5.*(pts[3].-centroid)
            fc=f(contract)
            if fc<min(fr,vals[3])
                pts[3]=contract; vals[3]=fc
            else
                for j in 2:3
                    pts[j]=pts[1].+0.5.*(pts[j].-pts[1]); vals[j]=f(pts[j])
                end
            end
        end
    end
    j=argmin(vals)
    return (pts[j], vals[j])
end
norm2(v)=sqrt(sum(abs2,v))

function binom_logp(c::Cell,p)
    return c.logcomb + c.k*log(p) + (c.n-c.k)*log1p(-p)
end
function bb_logp(c::Cell,p,phi)
    a=p*phi; b=(1-p)*phi
    return c.logcomb + logbeta(c.k+a,c.n-c.k+b) - logbeta(a,b)
end

struct Fit
    model::String
    p::Float64
    phi::Float64
    loglik::Float64
end
function fit_binomial(cells)
    total_k=sum(c.k*c.count for c in cells)
    total_n=sum(c.n*c.count for c in cells)
    total_n>0 || error("No opportunities")
    p=clamp_prob(total_k/total_n)
    return Fit("Binomial",p,Inf,sum(c.count*binom_logp(c,p) for c in cells))
end
function fit_betabinomial(cells)
    init=fit_binomial(cells)
    # Restarts are important because near-binomial data produce very flat likelihoods.
    bestval=Inf; bestx=[logit(init.p),log(5.0)]
    objective=function(x)
        (abs(x[1])>23 || x[2]<-14 || x[2]>13) && return 1e200
        p=clamp_prob(sigmoid(x[1])); phi=exp(x[2])
        ll=0.0
        # Precompute shared beta term once per likelihood evaluation.
        a=p*phi; b=(1-p)*phi; norm=logbeta(a,b)
        for c in cells
            ll+=c.count*(c.logcomb+logbeta(c.k+a,c.n-c.k+b)-norm)
        end
        return isfinite(ll) ? -ll : 1e200
    end
    for concentration in (0.05, 1.0, 20.0, 1000.0)
        x,val=nm_min(objective,[logit(init.p),log(concentration)])
        if val<bestval
            bestval=val; bestx=x
        end
    end
    return Fit("Beta-binomial",clamp_prob(sigmoid(bestx[1])),exp(bestx[2]),-bestval)
end

function logsumexp(v)
    m=maximum(v)
    return m+log(sum(exp(x-m) for x in v))
end

# Positive-support descriptive candidates, each normalized on k=1,...,N.
function log_descriptive_weights(model,N,param)
    if model=="Truncated geometric"
        lambda=param[1]
        return [-lambda*(k-1) for k in 1:N]
    elseif model=="Truncated power law"
        alpha=param[1]
        return [-alpha*log(k) for k in 1:N]
    elseif model=="Shifted truncated negative binomial"
        r=exp(param[1]); p=clamp_prob(sigmoid(param[2]))
        return [(
            loggamma_pos(k - 1 + r) - loggamma_pos(r) - loggamma_pos(k) +
            r * log(p) + (k - 1) * log1p(-p)
        ) for k in 1:N]
    end
    error("Unrecognised descriptive model")
end
function fit_descriptive(data,out)
    rows=NamedTuple[]; curves=NamedTuple[]
    for ds in DATA
        v=data[ds]; N=Int(v.N[1]); hist=zeros(Int,N)
        for k in v.K
            k>0 && (hist[Int(k)]+=1)
        end
        ns=sum(hist); empirical=hist./ns
        for model in ("Truncated geometric","Truncated power law","Shifted truncated negative binomial")
            loss=function(x)
                if any(y->!isfinite(y),x) || (model=="Shifted truncated negative binomial" &&
                    (x[1]<-14 || x[1]>12 || abs(x[2])>20))
                    return 1e200
                end
                w=log_descriptive_weights(model,N,x)
                z=logsumexp(w)
                ans=-sum(hist[k]*(w[k]-z) for k in 1:N if hist[k]>0)
                return isfinite(ans) ? ans : 1e200
            end
            if model=="Truncated geometric"
                theta,llneg=golden_max(lambda -> -loss([lambda]),0.0,16.0)
                pars=[theta]; nparam=1; nll=-llneg
            elseif model=="Truncated power law"
                theta,llneg=golden_max(alpha -> -loss([alpha]),0.001,12.0)
                pars=[theta]; nparam=1; nll=-llneg
            else
                mu=max(mean(Float64.(v.K[v.K .> 0]))-1,0.01)
                bestval=Inf; bestx=[0.0,0.0]
                for r0 in (0.2,1.0,5.0)
                    p0=clamp_prob(r0/(r0+mu))
                    x,val=nm_min(loss,[log(r0),logit(p0)])
                    if val<bestval
                        bestval=val; bestx=x
                    end
                end
                pars=bestx; nparam=2; nll=bestval
            end
            w=log_descriptive_weights(model,N,pars)
            q=exp.(w .- logsumexp(w))
            maxcdf=maximum(abs.(cumsum(q).-cumsum(empirical)))
            for k in 1:N
                push!(curves,(dataset=ds,model=model,K=k,expected_probability=q[k],
                              observed_probability=empirical[k]))
            end
            push!(rows,(dataset=ds,model=model,n_positive=ns,parameter_1=pars[1],
                        parameter_2=length(pars)==2 ? pars[2] : missing,
                        loglik=-nll,n_parameters=nparam,AIC=2*nparam+2*nll,
                        BIC=log(ns)*nparam+2*nll,KS_cdf=maxcdf,
                        note="Descriptive fit conditional on K>0; bounded support 1:N"))
        end
    end
    scores=DataFrame(rows)
    scores.delta_AIC=zeros(Float64,nrow(scores))
    for ds in DATA
        ix=findall(==(ds),scores.dataset)
        a=minimum(scores.AIC[ix]); scores.delta_AIC[ix]=scores.AIC[ix].-a
    end
    dir=joinpath(out,"02_distribution_fits")
    store(joinpath(dir,"descriptive_model_comparison.csv"),scores)
    store(joinpath(dir,"descriptive_fitted_pmf.csv"),DataFrame(curves))
    return scores,DataFrame(curves)
end

function by_n_counts(v)
    ans=Dict{Int,Int}()
    for n in v.n
        ni=Int(n); ans[ni]=get(ans,ni,0)+1
    end
    return ans
end

function pmf_fixed_n(n::Int, fit::Fit)
    # Log-space normalization prevents underflow for large site counts / small p.
    p=fit.p; q=zeros(Float64,n+1)
    if fit.model=="Binomial"
        logw=[logchoose(n,k)+k*log(p)+(n-k)*log1p(-p) for k in 0:n]
    else
        a=p*fit.phi; b=(1-p)*fit.phi; norm=logbeta(a,b)
        logw=[logchoose(n,k)+logbeta(k+a,n-k+b)-norm for k in 0:n]
    end
    z=logsumexp(logw)
    for k in 0:n
        q[k+1]=exp(logw[k+1]-z)
    end
    return q
end
function numerical_self_check()
    for fit in (Fit("Binomial",0.31,Inf,0.0),Fit("Beta-binomial",0.31,3.8,0.0))
        q=pmf_fixed_n(12,fit)
        @assert isapprox(sum(q),1.0;atol=1e-11)
        @assert isapprox(sum((k-1)*q[k] for k in eachindex(q)),12*fit.p;atol=1e-9)
    end
end
function expected_from_opportunities(v,fit)
    N=Int(v.N[1]); counts=by_n_counts(v); expected=zeros(Float64,N+1)
    for (n,c) in counts
        q=pmf_fixed_n(n,fit)
        @views expected[1:n+1].+=c.*q
    end
    return expected./nrow(v)
end
function empirical_probs(v)
    N=Int(v.N[1]); counts=zeros(Int,N+1)
    for k in v.K
        counts[Int(k)+1]+=1
    end
    return counts./nrow(v)
end
function js_divergence(p,q)
    m=(p.+q)./2
    return 0.5*sum(p[i]>0 ? p[i]*log(p[i]/m[i]) : 0.0 for i in eachindex(p))+
           0.5*sum(q[i]>0 ? q[i]*log(q[i]/m[i]) : 0.0 for i in eachindex(q))
end
function scoring_row(ds,context,fit,v)
    observed=empirical_probs(v)
    predicted=expected_from_opportunities(v,fit)
    pairs=nrow(v); N=Int(v.N[1])
    cells=cells_from(v)
    ll = if fit.model == "Binomial"
        sum(c.count * binom_logp(c, fit.p) for c in cells)
    else
        sum(c.count * bb_logp(c, fit.p, fit.phi) for c in cells)
    end
    return (dataset=ds,evaluation=context,model=fit.model,p=fit.p,
            concentration=fit.model=="Binomial" ? missing : fit.phi,
            n_pairs=pairs,mean_negative_loglik=-ll/pairs,
            JS_divergence=js_divergence(observed,predicted),
            total_variation=0.5*sum(abs.(observed.-predicted)),
            KS_CDF=maximum(abs.(cumsum(observed).-cumsum(predicted))),
            zero_error_pp=100*(predicted[1]-observed[1]),
            singleton_error_pp=100*(predicted[2]-observed[2]),
            singleton_positive_error_pp=100*(predicted[2]/max(1-predicted[1],1e-12)-
                                              observed[2]/max(1-observed[1],1e-12)),
            tail_ge_5_error_pp=N>=5 ? 100*(sum(predicted[6:end])-sum(observed[6:end])) : missing)
end
function fit_and_validate(data,out)
    cells=Dict(ds=>cells_from(data[ds]) for ds in DATA)
    pooled=merge_cells([cells[ds] for ds in DATA])
    pooled_bin=fit_binomial(pooled)
    pooled_bb=fit_betabinomial(pooled)
    fits=NamedTuple[]; scores=NamedTuple[]; prediction=NamedTuple[]
    local_fits=Dict{String,Tuple{Fit,Fit}}()
    loo_fits=Dict{String,Tuple{Fit,Fit}}()
    function record_fit(scope,heldout,fit,nobs)
        push!(fits,(scope=scope,dataset=heldout,model=fit.model,
                    p=fit.p,concentration=fit.model=="Binomial" ? missing : fit.phi,
                    loglik_training=fit.loglik,n_training_pairs=nobs,
                    AIC=2*(fit.model=="Binomial" ? 1 : 2)-2*fit.loglik))
    end
    record_fit("pooled_all","ALL",pooled_bin,sum(c.count for c in pooled))
    record_fit("pooled_all","ALL",pooled_bb,sum(c.count for c in pooled))
    for ds in DATA
        local_bin=fit_binomial(cells[ds]); local_bb=fit_betabinomial(cells[ds])
        local_fits[ds]=(local_bin,local_bb)
        record_fit("dataset_only",ds,local_bin,nrow(data[ds]))
        record_fit("dataset_only",ds,local_bb,nrow(data[ds]))
        train=merge_cells([cells[other] for other in DATA if other != ds])
        loo_bin=fit_binomial(train); loo_bb=fit_betabinomial(train)
        loo_fits[ds]=(loo_bin,loo_bb)
        record_fit("leave_one_dataset_out",ds,loo_bin,sum(c.count for c in train))
        record_fit("leave_one_dataset_out",ds,loo_bb,sum(c.count for c in train))
        for (context, pair) in (("dataset_fitted",(local_bin,local_bb)),
                                ("pooled_fitted",(pooled_bin,pooled_bb)),
                                ("leave_one_dataset_out",(loo_bin,loo_bb)))
            for fit in pair
                push!(scores,scoring_row(ds,context,fit,data[ds]))
                emp=empirical_probs(data[ds]); expected=expected_from_opportunities(data[ds],fit)
                for k in 0:Int(data[ds].N[1])
                    push!(prediction,(dataset=ds,evaluation=context,model=fit.model,K=k,
                                      observed=emp[k+1],expected=expected[k+1],
                                      observed_ccdf=sum(emp[k+1:end]),expected_ccdf=sum(expected[k+1:end])))
                end
            end
        end
        println("  LOO complete: ",ds," (BB p=",round(loo_bb.p,digits=4),
                ", phi=",round(loo_bb.phi,digits=3),")")
    end
    dir=joinpath(out,"03_opportunity_models")
    store(joinpath(dir,"model_parameters.csv"),DataFrame(fits))
    store(joinpath(dir,"model_predictions_pmf_ccdf.csv"),DataFrame(prediction))
    scoresdf=DataFrame(scores)
    store(joinpath(dir,"model_predictive_scores.csv"),scoresdf)
    loodf=scoresdf[scoresdf.evaluation.=="leave_one_dataset_out",:]
    store(joinpath(out,"04_cross_dataset_generality","leave_one_dataset_out_scores.csv"),loodf)
    pooledsummary=combine(groupby(loodf,:model),
                          :mean_negative_loglik=>mean=>:equal_dataset_mean_pair_NLL,
                          :JS_divergence=>mean=>:mean_dataset_JS,
                          :total_variation=>mean=>:mean_dataset_TV,
                          :KS_CDF=>mean=>:mean_dataset_KS,
                          :zero_error_pp=>(x->mean(abs.(x)))=>:mean_absolute_zero_error_pp,
                          :singleton_error_pp=>(x->mean(abs.(x)))=>:mean_absolute_singleton_error_pp,
                          :singleton_positive_error_pp=>(x->mean(abs.(x)))=>:mean_absolute_positive_singleton_error_pp)
    store(joinpath(out,"04_cross_dataset_generality","loo_comparison_summary.csv"),pooledsummary)
    # Quantify pooling-versus-heterogeneity; parameter counts account for all 10 datasets.
    comparisons=NamedTuple[]
    for (name,shared,k_per_ds) in (("Binomial",pooled_bin,1),("Beta-binomial",pooled_bb,2))
        local_ll=sum(name=="Binomial" ? local_fits[d][1].loglik : local_fits[d][2].loglik for d in DATA)
        common_aic=2*k_per_ds-2*shared.loglik
        separate_aic=2*k_per_ds*length(DATA)-2*local_ll
        push!(comparisons,(model=name,shared_loglik=shared.loglik,separate_loglik=local_ll,
                           n_datasets=length(DATA),shared_parameters=k_per_ds,
                           separate_parameters=k_per_ds*length(DATA),
                           shared_AIC=common_aic,separate_AIC=separate_aic,
                           AIC_separate_minus_shared=separate_aic-common_aic,
                           note="In-sample descriptive AIC; pair/site dependence can inflate evidence"))
    end
    store(joinpath(out,"04_cross_dataset_generality","pooled_vs_dataset_specific_fits.csv"),DataFrame(comparisons))
    return local_fits,loo_fits,scoresdf,DataFrame(prediction)
end

function simulate_envelopes(data, loo_fits, out, replicates, rng)
    replicates==0 && return DataFrame()
    envelopes=NamedTuple[]
    for ds in DATA
        v=data[ds]; N=Int(v.N[1]); byn=by_n_counts(v); npairs=nrow(v)
        for fit in loo_fits[ds]
            dist=Dict(n=>cumsum(pmf_fixed_n(n,fit)) for n in keys(byn))
            for cumulative in values(dist)
                cumulative[end]=1.0  # protect categorical sampling against rounding
            end
            simulations=zeros(Float64,N+1,replicates)
            for rep in 1:replicates
                counts=zeros(Int,N+1)
                for (n,c) in byn
                    cumulative=dist[n]
                    for _ in 1:c
                        k=searchsortedfirst(cumulative,rand(rng))-1
                        counts[min(k,N)+1]+=1
                    end
                end
                simulations[:,rep]=counts./npairs
            end
            observed=empirical_probs(v)
            pred=expected_from_opportunities(v,fit)
            # CCDF is P(K>=k); for k=0 the full distribution equals 1.
            tails=reverse(cumsum(reverse(simulations,dims=1),dims=1),dims=1)
            for k in 0:N
                push!(envelopes,(dataset=ds,model=fit.model,K=k,
                      observed=observed[k+1],expected=pred[k+1],
                      pmf_lower=quantile(view(simulations,k+1,:),0.025),
                      pmf_upper=quantile(view(simulations,k+1,:),0.975),
                      ccdf_observed=sum(observed[k+1:end]),
                      ccdf_expected=sum(pred[k+1:end]),
                      ccdf_lower=quantile(view(tails,k+1,:),0.025),
                      ccdf_upper=quantile(view(tails,k+1,:),0.975),
                      replicates=replicates))
            end
        end
    end
    store(joinpath(out,"03_opportunity_models","loo_parametric_envelopes.csv"),DataFrame(envelopes))
    return DataFrame(envelopes)
end

function do_site_subsampling(data,input,out,reps,rng)
    reps==0 && return (DataFrame(),DataFrame())
    target=minimum(Int(data[ds].N[1]) for ds in DATA)
    if target<3
        @warn "Site subsampling skipped: smallest dataset has fewer than 3 sites" target
        return (DataFrame(),DataFrame())
    end
    summaries=NamedTuple[]; curves=NamedTuple[]
    for ds in DATA
        occ=CSV.read(joinpath(input,"inputs",ds*"_occupancy.csv"),DataFrame;types=String)
        sites=sort(unique(occ.site))
        co=CSV.read(joinpath(input,"inputs",ds*"_cooccurrences.csv"),DataFrame;types=String)
        ints=CSV.read(joinpath(input,"inputs",ds*"_interactions.csv"),DataFrame;types=String)
        co_by_site=Dict(s=>Tuple{String,String}[] for s in sites)
        int_by_site=Dict(s=>Tuple{String,String}[] for s in sites)
        for r in eachrow(co)
            push!(co_by_site[r.site],(r.consumer,r.resource))
        end
        for r in eachrow(ints)
            push!(int_by_site[r.site],(r.consumer,r.resource))
        end
        for rep in 1:reps
            chosen=sites[randperm(rng,length(sites))[1:target]]
            n=Dict{Tuple{String,String},Int}(); k=Dict{Tuple{String,String},Int}()
            for s in chosen
                for pair in co_by_site[s]
                    n[pair]=get(n,pair,0)+1
                end
                for pair in int_by_site[s]
                    k[pair]=get(k,pair,0)+1
                end
            end
            all_count=length(n); positive=count(x->x>0,values(k)); singletons=count(x->x==1,values(k))
            all_count>0 || continue
            hist=zeros(Int,target+1)
            for pair in keys(n)
                hist[get(k,pair,0)+1]+=1
            end
            push!(summaries,(dataset=ds,replicate=rep,sites_sampled=target,
                             candidate_pairs=all_count,positive_pairs=positive,
                             zero_fraction=hist[1]/all_count,
                             singleton_fraction_positive=positive>0 ? singletons/positive : missing,
                             mean_K_positive=positive>0 ? sum(values(k))/positive : missing))
            for t in 0:target
                push!(curves,(dataset=ds,replicate=rep,sites_sampled=target,K=t,
                              probability=hist[t+1]/all_count,
                              ccdf=sum(hist[t+1:end])/all_count))
            end
        end
    end
    dir=joinpath(out,"05_site_subsampling")
    store(joinpath(dir,"site_subsampling_summary.csv"),DataFrame(summaries))
    store(joinpath(dir,"site_subsampling_curves.csv"),DataFrame(curves))
    return (DataFrame(summaries),DataFrame(curves))
end

function configure_theme()
    set_theme!(Theme(font="Arial",fontsize=18,
        Axis=(topspinevisible=false,rightspinevisible=false,
              xgridvisible=false,ygridvisible=false,
              xticklabelsize=13,yticklabelsize=13)))
end
function panel_data(df,ds;positive_only=nothing)
    idx=df.dataset .== ds
    if positive_only !== nothing
        idx .&= df.positive_only .== positive_only
    end
    return df[idx,:]
end

function draw_empirical(pmf,opp,out)
    dir=joinpath(out,"01_empirical_distributions")
    for (metric,label,name) in ((:probability,"P(K=k | K>0)","01_positive_support_pmf.png"),
                                (:ccdf,"P(K>=k | K>0)","02_positive_support_ccdf.png"))
        fig=Figure(size=(1900,800))
        for (i,ds) in enumerate(DATA)
            t=panel_data(pmf,ds;positive_only=true)
            t=t[(t.K .>= 1).&(t[!,metric].>0),:]
            ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
                    xlabel="Local support K",ylabel=mod1(i,5)==1 ? label : "",
                    xscale=log10,yscale=log10)
            scatterlines!(ax,t.K,t[!,metric],color=COLORS[i],markersize=5,linewidth=1.6)
        end
        save(joinpath(dir,name),fig;px_per_unit=2)
    end
    fig=Figure(size=(1750,620))
    ax1=Axis(fig[1,1],title="Positive support PMF",xlabel="K",ylabel="P(K=k | K>0)",xscale=log10,yscale=log10)
    ax2=Axis(fig[1,2],title="Positive support CCDF",xlabel="K",ylabel="P(K>=k | K>0)",xscale=log10,yscale=log10)
    ax3=Axis(fig[1,3],title="Opportunity CCDF",xlabel="n",ylabel="P(n>=x)",xscale=log10,yscale=log10)
    for (i,ds) in enumerate(DATA)
        t=panel_data(pmf,ds;positive_only=true)
        t=t[(t.K.>=1).&(t.probability.>0),:]
        scatterlines!(ax1,t.K,t.probability,color=(COLORS[i],0.8),markersize=3,linewidth=1.5)
        scatterlines!(ax2,t.K,t.ccdf,color=(COLORS[i],0.8),markersize=3,linewidth=1.5)
        o=opp[opp.dataset.==ds,:]
        o=o[(o.n.>=1).&(o.ccdf.>0),:]
        lines!(ax3,o.n,o.ccdf,color=(COLORS[i],0.8),linewidth=2)
    end
    Legend(fig[2,1:3],[LineElement(color=c,linewidth=3) for c in COLORS],
           replace.(DATA,"_"=>" ");orientation=:horizontal,nbanks=2,framevisible=false,labelsize=15)
    save(joinpath(dir,"03_overlay_support_and_opportunities.png"),fig;px_per_unit=2)
    # Full candidate-pair distribution including zeros (linear axes preserve zero spike).
    fig0=Figure(size=(1850,750))
    for (i,ds) in enumerate(DATA)
        t=panel_data(pmf,ds;positive_only=false)
        t=t[t.count.>0,:]
        ax=Axis(fig0[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
                xlabel="Local support K",ylabel=mod1(i,5)==1 ? "P(K=k)" : "")
        scatterlines!(ax,t.K,t.probability;color=COLORS[i],linewidth=1.5,markersize=5)
    end
    save(joinpath(dir,"04_full_pmf_including_zeros.png"),fig0;px_per_unit=2)
end

function draw_conditional_heatmap(out)
    raw=CSV.read(joinpath(out,"01_empirical_distributions","conditional_fraction_by_n_bin.csv"),DataFrame)
    fig=Figure(size=(1900,870))
    nlabels=("1", "2-3", "4-7", "8-15", "16-31", "32-63", "64+")
    flabels=("0", "(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")
    for (i,ds) in enumerate(DATA)
        ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
            xlabel="Opportunity n (binned)",ylabel=mod1(i,5)==1 ? "K/n category" : "",
            xticks=(1:7,collect(nlabels)),yticks=(1:5,collect(flabels)),
            xticklabelrotation=pi/3,xticklabelsize=10,yticklabelsize=10)
        matrix=fill(NaN,7,5)
        t=raw[(raw.dataset.==ds).&(.!raw.positive_only),:]
        for row in eachrow(t)
            x=findfirst(==(row.n_bin),nlabels)
            y=findfirst(==(row.fraction_bin),flabels)
            (x===nothing || y===nothing) && continue
            matrix[x,y]=row.probability
        end
        heatmap!(ax,1:7,1:5,matrix,colormap=:viridis,colorrange=(0,1))
    end
    Colorbar(fig[1:2,6],colormap=:viridis,colorrange=(0,1),label="Conditional probability")
    save(joinpath(out,"01_empirical_distributions","05_conditional_realisation_heatmaps.png"),fig;px_per_unit=2)
end

function draw_descriptive_fits(curves,out)
    fig=Figure(size=(1900,850))
    models=("Truncated geometric","Truncated power law","Shifted truncated negative binomial")
    cols=("#0072B2","#D55E00","#009E73")
    for (i,ds) in enumerate(DATA)
        ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
                xlabel="Local support K",ylabel=mod1(i,5)==1 ? "Positive-support PMF" : "",
                xscale=log10,yscale=log10)
        emp=curves[(curves.dataset.==ds).&(curves.model.==models[1]),:]
        emp=emp[emp.observed_probability.>0,:]
        scatter!(ax,emp.K,emp.observed_probability,color=:black,markersize=7)
        for (m,col) in zip(models,cols)
            t=curves[(curves.dataset.==ds).&(curves.model.==m),:]
            t=t[t.expected_probability.>1e-12,:]
            lines!(ax,t.K,t.expected_probability;color=col,linewidth=2)
        end
    end
    items=Any[MarkerElement(color=:black,marker=:circle,markersize=11)]
    append!(items,[LineElement(color=c,linewidth=3) for c in cols])
    Legend(fig[3,1:5],items,["Observed";collect(models)],orientation=:horizontal,
           nbanks=1,framevisible=false,labelsize=17)
    save(joinpath(out,"02_distribution_fits","01_descriptive_candidate_fits.png"),fig;px_per_unit=2)
end

function draw_model_predictions(prediction,out)
    dir=joinpath(out,"03_opportunity_models")
    fig=Figure(size=(1900,850))
    for (i,ds) in enumerate(DATA)
        ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
                xlabel="Local support K (showing 0-40)",
                ylabel=mod1(i,5)==1 ? "Probability" : "")
        for (model,col) in (("Binomial","#D55E00"),("Beta-binomial","#0072B2"))
            t=prediction[(prediction.dataset.==ds).&
                         (prediction.evaluation.=="leave_one_dataset_out").&
                         (prediction.model.==model),:]
            t=t[(t.expected.>0).&(t.K.<=min(40,maximum(t.K))),:]
            lines!(ax,t.K,t.expected;color=col,linewidth=2.3)
        end
        t=prediction[(prediction.dataset.==ds).&
                     (prediction.evaluation.=="leave_one_dataset_out").&
                     (prediction.model.=="Binomial"),:]
        t=t[(t.observed.>0).&(t.K.<=40),:]
        scatter!(ax,t.K,t.observed;color=:black,markersize=6)
    end
    legenditems=[MarkerElement(color=:black,marker=:circle,markersize=12),
                 LineElement(color="#D55E00",linewidth=3),LineElement(color="#0072B2",linewidth=3)]
    Legend(fig[3,1:5],legenditems,["Observed","LOO binomial","LOO beta-binomial"],
           orientation=:horizontal,framevisible=false)
    save(joinpath(dir,"01_loo_opportunity_controlled_pmf.png"),fig;px_per_unit=2)
end

function draw_envelopes(envelopes,out)
    nrow(envelopes)==0 && return
    dir=joinpath(out,"03_opportunity_models")
    fig=Figure(size=(1850,850))
    for (i,ds) in enumerate(DATA)
        t=envelopes[(envelopes.dataset.==ds).&(envelopes.model.=="Beta-binomial").&(envelopes.K.>=1),:]
        ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),
                xlabel="Local support K",ylabel=mod1(i,5)==1 ? "P(K>=k)" : "",
                xscale=log10,yscale=log10)
        # A log-axis cannot render zero CI bounds: clip only the plotted lower envelope.
        lower=max.(t.ccdf_lower,1e-7)
        upper=max.(t.ccdf_upper,lower)
        band!(ax,t.K,lower,upper,color=("#0072B2",0.15))
        lines!(ax,t.K,max.(t.ccdf_expected,1e-7),color="#0072B2",linewidth=2.5)
        scatter!(ax,t.K[t.ccdf_observed.>0],t.ccdf_observed[t.ccdf_observed.>0],
                 color=:black,markersize=5)
    end
    save(joinpath(dir,"02_loo_beta_binomial_ccdf_envelopes.png"),fig;px_per_unit=2)
end
function draw_loo_scores(scores,out)
    t=scores[scores.evaluation.=="leave_one_dataset_out",:]
    fig=Figure(size=(1250,570))
    ax=Axis(fig[1,1],title="Out-of-dataset predictive performance",
            ylabel="Mean negative log-likelihood per pair",xlabel="Held-out dataset",
            xticks=(1:10,replace.(DATA,"_"=>" ")),xticklabelrotation=pi/4)
    for (model,col,offset) in (("Binomial","#D55E00",-0.11),("Beta-binomial","#0072B2",0.11))
        vals=[only(t.mean_negative_loglik[(t.dataset.==ds).&(t.model.==model)]) for ds in DATA]
        scatter!(ax,collect(1:10).+offset,vals,color=col,markersize=14)
    end
    items=[MarkerElement(color="#D55E00",marker=:circle,markersize=12),
           MarkerElement(color="#0072B2",marker=:circle,markersize=12)]
    Legend(fig[2,1],items,["Binomial","Beta-binomial"],orientation=:horizontal,framevisible=false)
    save(joinpath(out,"04_cross_dataset_generality","01_loo_predictive_nll.png"),fig;px_per_unit=2)
end
function draw_subsampling(summaries,out)
    nrow(summaries)==0 && return
    fig=Figure(size=(1300,590))
    ax=Axis(fig[1,1],title="Regional link singleton fraction after equal-site subsampling",
            xlabel="Dataset",ylabel="Fraction of positive pairs with K=1",
            xticks=(1:10,replace.(DATA,"_"=>" ")),xticklabelrotation=pi/4)
    for (i,ds) in enumerate(DATA)
        vals=collect(skipmissing(summaries.singleton_fraction_positive[summaries.dataset.==ds]))
        isempty(vals) && continue
        scatter!(ax,fill(i,length(vals)),vals,color=(COLORS[i],0.22),markersize=6)
        scatter!(ax,[i],[mean(vals)],color=COLORS[i],markersize=17)
    end
    save(joinpath(out,"05_site_subsampling","01_singletons_after_equal_site_subsampling.png"),fig;px_per_unit=2)
end

function write_notes(out,cfg,minimum_sites)
    notes="""LOCAL-SUPPORT DISTRIBUTION ANALYSIS\n\n"""*
    """Input: outputs/support_prediction_zero/$(cfg.scenario)/, validated against site-level tables; stored SHA-256 file hashes. Source SHA checks are enforced when input_versions.csv is present.\n"""*
    """Primary inference: all co-occurring pairs, K>=0. Fit binomial and beta-binomial K|n using exact observed opportunity counts.\n"""*
    """Beta-binomial p is the across-pair mean realisation probability; phi is concentration (alpha=p*phi, beta=(1-p)*phi).\n"""*
    """LOO: parameters trained on nine whole datasets and applied without recalibration to the tenth; target n is supplied, target K is never used to fit.\n"""*
    """Likelihood uses independent pair-level contributions. Shared sites/species can create dependence: predictive scores are useful, but naive standard errors and AIC certainty would be optimistic.\n"""*
    """Empirical PMFs include zeros or condition on K>0. Descriptive families (finite geometric, power law, shifted negative binomial) condition on K>0 and truncate at 1:N.\n"""*
    """Power-law fits are descriptive candidates only, NOT proof of power-law behaviour or a formal Clauset-style tail test.\n"""*
    """Descriptive parameter_1: geometric lambda, power-law alpha, or negative-binomial log(shape). Descriptive parameter_2: negative-binomial logit(success probability).\n"""*
    """LOO parametric envelopes hold n, model parameters and candidate pairs fixed. They exclude model-parameter uncertainty and site/species dependence.\n"""*
    """Site subsampling uses $(cfg.subsamples) replicates at min dataset site count=$minimum_sites (unless disabled). It does not standardise spatial grain or study extent.\n"""*
    """Zero support means unobserved interaction among an inferred co-occurring pair, not necessarily a confirmed biological absence.\n"""*
    """Binomial and beta-binomial models abstract from network degree structure, species identities, spatial correlations, and detection differences.\n"""*
    """Simulations use seed=$(cfg.seed) and $(cfg.replicates) replicates per dataset per LOO model.\n"""*
    """All CSV and PNG outputs are generated inside the named scenario output folder.\n"""
    open(joinpath(out,"README_results.txt"),"w") do io
        write(io,notes)
    end
end

function main()
    cfg=parse_args(ARGS)
    input=joinpath(cfg.root,"outputs","support_prediction_zero",cfg.scenario)
    out=joinpath(cfg.root,"outputs","statistical description of local support",
                 "local_support_distributions",cfg.scenario)
    input_dir=joinpath(input,"inputs")
    isdir(input_dir) || error("Missing scenario inputs at $input_dir. Run prepare_zero_presence preprocessing first.")
    for d in ("01_empirical_distributions","02_distribution_fits","03_opportunity_models",
              "04_cross_dataset_generality","05_site_subsampling")
        mkpath(joinpath(out,d))
    end
    if cfg.replicates==0
        for f in (joinpath(out,"03_opportunity_models","loo_parametric_envelopes.csv"),
                  joinpath(out,"03_opportunity_models","02_loo_beta_binomial_ccdf_envelopes.png"))
            isfile(f) && rm(f)
        end
    end
    if cfg.subsamples==0
        for f in (joinpath(out,"05_site_subsampling","site_subsampling_summary.csv"),
                  joinpath(out,"05_site_subsampling","site_subsampling_curves.csv"),
                  joinpath(out,"05_site_subsampling","01_singletons_after_equal_site_subsampling.png"))
            isfile(f) && rm(f)
        end
    end
    rng=MersenneTwister(cfg.seed)
    numerical_self_check()
    println("Local support distributions: scenario=",cfg.scenario,"; destination=",out)
    pairs,pairs_existing=load_or_build_pairs(input,out)
    verify_data(input,pairs)
    data=Dict(ds=>pairs[pairs.dataset.==ds,:] for ds in DATA)
    manifest=NamedTuple[]
    paths=String[]
    pairs_existing && push!(paths,joinpath(input,"pairs.csv"))
    isfile(joinpath(input,"input_versions.csv")) && push!(paths,joinpath(input,"input_versions.csv"))
    for ds in DATA,kind in ("interactions","cooccurrences","occupancy")
        push!(paths,joinpath(input,"inputs",ds*"_"*kind*".csv"))
    end
    for path in paths
        push!(manifest,(file=relpath(path,cfg.root),sha256=file_sha(path)))
    end
    push!(manifest,(file=relpath(@__FILE__,cfg.root),sha256=file_sha(@__FILE__)))
    store(joinpath(out,"input_and_code_checksums.csv"),DataFrame(manifest))
    pmf,opportunities,summary=basic_distributions(data,out)
    println("  descriptive distributions complete")
    descriptive_scores,curves=fit_descriptive(data,out)
    println("  descriptive candidate fits complete")
    localfits,loofits,predictive_scores,prediction=fit_and_validate(data,out)
    envelopes=simulate_envelopes(data,loofits,out,cfg.replicates,rng)
    subsampling_summary,subsampling_curves=do_site_subsampling(data,input,out,cfg.subsamples,rng)
    configure_theme()
    draw_empirical(pmf,opportunities,out)
    draw_conditional_heatmap(out)
    draw_descriptive_fits(curves,out)
    draw_model_predictions(prediction,out)
    draw_envelopes(envelopes,out)
    draw_loo_scores(predictive_scores,out)
    draw_subsampling(subsampling_summary,out)
    write_notes(out,cfg,minimum(Int(data[ds].N[1]) for ds in DATA))
    println("COMPLETE. Outputs written to: ",out)
    println("  LOO summary: ",joinpath(out,"04_cross_dataset_generality","loo_comparison_summary.csv"))
end

main()
end # module
