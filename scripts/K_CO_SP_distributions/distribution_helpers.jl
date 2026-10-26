function normprob(logw)
    mx=maximum(logw);w=exp.(logw.-mx);z=sum(w)
    w/z,mx+log(z)
end
function onefit(x,target)
    p,_=normprob(zeros(length(x)))
    target>=dot(p,x)-1e-13 && return 0.
    target<=minimum(x)+1e-13 && return 1000.
    lo=0.;hi=1.
    while dot(first(normprob(-hi.*x)),x)>target && hi<1e5;hi*=2;end
    for _ in 1:48
        mid=(lo+hi)/2
        if dot(first(normprob(-mid.*x)),x)>target;lo=mid;else;hi=mid;end
    end
    (lo+hi)/2
end
function objective(theta,X,target)
    p,z=normprob(-X*theta);mu=transpose(X)*p;g=target.-mu
    H=transpose(X)*(X.*p).-mu*transpose(mu)
    dot(theta,target)+z,g,H,p
end
function optimise(X,target,start,positive)
    theta=copy(start);best=Inf
    for iter in 1:150
        value,g,H,p=objective(theta,X,target)
        free=[j for j in eachindex(theta) if !positive[j] || theta[j]>1e-10 || g[j]<0]
        isempty(free) && break
        maximum(abs.(g[free]))<1e-9 && break
        delta=zeros(length(theta));delta[free]=-(H[free,free]+1e-12I)\g[free]
        dot(delta,g)>=0 && (delta[free]=-g[free])
        step=1.;accepted=false
        for _ in 1:45
            proposal=theta.+step.*delta
            for j in eachindex(theta);positive[j] && (proposal[j]=max(0.,proposal[j]));end
            nv=objective(proposal,X,target)[1]
            if nv<=value+1e-4dot(g,proposal.-theta)+1e-14
                theta=proposal;best=nv;accepted=true;break
            end
            step/=2
        end
        !accepted && break
    end
    val,g,_,p=objective(theta,X,target)
    residual=maximum([positive[j] && theta[j]<1e-8 ? max(0.,-g[j]) : abs(g[j]) for j in eachindex(theta)])
    (theta=theta,value=val,p=p,kkt=residual)
end
function fit(counts,model,xmin=1)
    N=length(counts);k=Float64.(xmin:N);n=sum(counts[xmin:end]);obs=counts[xmin:end]./n
    alpha=NaN;lambda=NaN;mu=NaN;sigma=NaN;s=NaN;r=NaN;beta=NaN;residual=0.;boundary=false
    if model=="Exponential"
        lambda=onefit(k,dot(obs,k));p,z=normprob(-lambda.*k);d=1;boundary=lambda==0
    elseif model=="Power_law"
        alpha=onefit(log.(k),dot(obs,log.(k)));p,z=normprob(-alpha.*log.(k));d=1;boundary=alpha==0
    elseif model in ["Cutoff_power_law","Discrete_lognormal"]
        lk=log.(k);alpha0=onefit(lk,dot(obs,lk))
        X=model=="Cutoff_power_law" ? hcat(lk,k) : hcat(lk,lk.^2)
        scale=[max(std(X[:,j]),1e-3) for j in 1:2];Xs=X./transpose(scale);target=transpose(Xs)*obs
        positive=model=="Cutoff_power_law" ? [true,true] : [false,true]
        starts=[scale.*[alpha0,0.],scale.*[alpha0,.1]]
        if model=="Cutoff_power_law";push!(starts,scale.*[0.,onefit(k,dot(obs,k))]);else;push!(starts,scale.*[-1.,1.]);end
        candidates=[optimise(Xs,target,start,positive) for start in starts];best=candidates[argmin([x.value for x in candidates])]
        theta=best.theta./scale;p=best.p;residual=best.kkt;alpha=theta[1];d=2
        if model=="Cutoff_power_law";lambda=theta[2];boundary=alpha<1e-8 || lambda<1e-8
        else;beta=theta[2];boundary=beta<1e-8;sigma=boundary ? Inf : sqrt(1/(2beta));mu=boundary ? NaN : (1-alpha)/(2beta);end
    else
        @assert xmin==1
        s=counts[1]/n;tail=counts[2:end];nt=sum(tail)
        lambda=nt==0 ? 1000. : onefit(Float64.(2:N),dot(tail./nt,Float64.(2:N)))
        pt,_=normprob(-lambda.*Float64.(2:N));p=vcat(s,(1-s).*pt);r=exp(-lambda);d=2;boundary=s==1 || lambda==0
    end
    @assert abs(sum(p)-1)<1e-10 && minimum(p)>=0
    ll=sum(counts[xmin+i-1]*log(p[i]) for i in eachindex(p) if counts[xmin+i-1]>0)
    D=maximum(abs.(cumsum(obs).-cumsum(p)));aic=2d-2ll;aicc=n>d+1 ? aic+2d*(d+1)/(n-d-1) : Inf
    (p=p,ll=ll,aic=aic,aicc=aicc,D=D,parameters=d,alpha=alpha,lambda=lambda,mu=mu,sigma=sigma,beta=beta,s=s,r=r,kkt=residual,boundary=boundary)
end
function samplecounts(rng,p,n)
    cdf=cumsum(p);cdf[end]=1.;counts=zeros(Int,length(p))
    for _ in 1:n;counts[searchsortedfirst(cdf,rand(rng))]+=1;end
    counts
end
function bootstrap(counts,model,f,rng)
    vals=zeros(B);maxkkt=0.
    for b in 1:B
        sim=samplecounts(rng,f.p,sum(counts));sf=fit(sim,model);vals[b]=sf.D;maxkkt=max(maxkkt,sf.kkt)
    end
    p=(1+count(>=(f.D-1e-12),vals))/(B+1)
    (p=p,mcse=sqrt(p*(1-p)/(B+1)),maxkkt=maxkkt)
end
