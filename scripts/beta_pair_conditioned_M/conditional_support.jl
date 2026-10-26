function gapmatrix(N)
    H=zeros(N+1,N)
    for m in 0:N,k in 1:N
        H[m+1,k]=m/N-(m<k ? 0. : prod((m-j)/(N-j) for j in 0:k-1))
    end
    H
end
function pairweights(C,a,h)
    [begin
        qs=logistic.([a-h,a+h])
        [sum(log((c-j+1)/j) for j in 1:k)+lse([k*log(q)+(c-k)*log1p(-q)-logz(c*log1p(-q)) for q in qs])-log(2) for k in 1:c]
    end for c in C]
end
function dp(logw,M)
    C=length.(logw);L=length(C)
    low=-60.;high=60.
    for _ in 1:65
        t=(low+high)/2
        mu=sum(begin v=w.+t.*collect(1:length(w));p=exp.(v.-lse(v));sum(p.*collect(1:length(w))) end for w in logw)
        if mu>M;high=t;else;low=t;end
    end
    tilt=(low+high)/2
    P=[exp.(w.+tilt.*collect(1:length(w)).-lse(w.+tilt.*collect(1:length(w)))) for w in logw]
    D=zeros(L+1,M+1);D[1,1]=1.;cap=0
    for i in 1:L
        cap+=C[i]
        for total in i:min(cap,M),k in 1:min(C[i],total-i+1)
            D[i+1,total+1]+=P[i][k]*D[i,total-k+1]
        end
        D[i+1,:]./=maximum(@view D[i+1,:])
    end
    @assert D[end,M+1]>0
    D,P
end
function drawK(rng,M,D,P)
    K=zeros(Int,length(P));left=M
    for i in length(P):-1:1
        w=[P[i][k]*D[i,left-k+1] for k in 1:min(length(P[i]),left-i+1)]
        u=rand(rng)*sum(w);v=0.;chosen=length(w)
        for k in eachindex(w);v+=w[k];if u<=v;chosen=k;break;end;end
        K[i]=chosen;left-=chosen
    end
    @assert left==0
    K
end
function exactfreq(M,D,P,N)
    state=zeros(M+1);state[M+1]=1.;freq=zeros(N)
    for i in length(P):-1:1
        prior=zeros(M+1)
        for total in i:M
            state[total+1]<1e-300 && continue
            upto=min(length(P[i]),total-i+1)
            den=sum(P[i][k]*D[i,total-k+1] for k in 1:upto)
            den==0 && continue
            for k in 1:upto
                mass=state[total+1]*P[i][k]*D[i,total-k+1]/den
                freq[k]+=mass;prior[total-k+1]+=mass
            end
        end
        state=prior
    end
    @assert isapprox(sum(freq),length(P);atol=1e-7)
    @assert isapprox(sum(freq.*collect(1:N)),M;atol=1e-7)
    freq/length(P)
end

