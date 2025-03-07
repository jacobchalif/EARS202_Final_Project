using Pkg
Pkg.activate(".")
Pkg.add(["StatGeochemBase","NaNStatistics","Distributions","Plots","DataFrames","DelimitedFiles"])
using StatGeochemBase, NaNStatistics, Distributions, Plots, DataFrames, DelimitedFiles

# Log likelihood functions for use within `align`
function yll(ymi, yn, y_sigma)
    ll = 0.0
    dist = Normal(0, y_sigma)
    @inbounds for i in eachindex(ymi, yn)
        ll += logpdf(dist, ymi[i]-yn[i])
    end
    return ll
end
function dll(d, x_sigma)
    ll = 0.0
    dist = Normal(0, x_sigma)
    @inbounds for i in eachindex(d)
        ll += logpdf(dist, d[i])
    end
    return ll
end

"""
```julia
align(xc, yc, xn, yn; 
    x_sigma::Number, 
    y_sigma::Number, 
    nsteps=10^5, 
    burnin=10^4
)
```
Align noisy signal `xn`,`yn` to central signal `xc`,`yc` by stretching `xn`,
using a Markov chain Monte Carlo approach given contstant uncertainties 
`x_sigma`, `y_sigma`, collecting `nsteps` after `burnin`. 

Returns results `ddist` (results), `lldist` (log likelihoods), `proposaldist` (proposal distribution sigma), and `acceptancedist` (boolean)
Each column of `ddist` represents a vector of perturbations `d` to `xn` drawn from the
stationary distribution of the markov chain, which can be applied as corrections to `xn`.

### Examples
```julia
julia> ddist, lldist, proposaldist, acceptancedist = align(xc, yc, xn, yn; x_sigma=1, y_sigma=0.1, nsteps=10^5, burnin=10^4);

julia> d_mean = nanmean(ddist, dim=2)
51-element Vector{Float64}:
 -0.3448032001182671
 -0.33213628128827566
  ⋮
  0.5484380519749825
```
"""
function align(xc, yc, xn, yn; x_sigma::Number, y_sigma::Number, nsteps=10^5, burnin=10^4)
    @assert issorted(xc) "Central x record must be sorted in increasing order"
    @assert issorted(xn) "Initial xn values must be sorted in increasing order"

    # Allocate variables for result distributions
    ddist = zeros(length(xn), nsteps)
    lldist = zeros(nsteps+burnin)
    proposaldist = zeros(nsteps)
    acceptancedist = falses(nsteps)

    # Allocate variables used in inversion
    d = zeros(size(xn))
    dₚ = zeros(size(xn))
    x1ₚ = copy(xn)
    ymi = similar(yn)
    jump_sigma = x_sigma/10


    # Likelihood of initial proposal
    linterp1!(ymi, xc, yc, x1ₚ)
    ll = llₚ = if issorted(x1ₚ)
        linterp1!(ymi, xc, yc, x1ₚ)
        yll(ymi, yn, y_sigma) + dll(dₚ, x_sigma)
    else
        -Inf
    end
    
    for i in 1:burnin
        # Reset proposals
        copyto!(dₚ, d)

        # Perturb displacement in one region using a Gaussian kernel
        r = rand(1:length(d))
        δ = randn()*jump_sigma
        bandwidth = rand()*x_sigma
        kernel = Normal(xn[r]+d[r], bandwidth)
        for i in eachindex(dₚ)
            dₚ[i] += δ * pdf(kernel, xn[i]+d[i]) * bandwidth
        end

        # Apply perturbation
        @. x1ₚ = xn + dₚ

        # Calculate log likelihood
        llₚ = if issorted(x1ₚ)
            linterp1!(ymi, xc, yc, x1ₚ)
            yll(ymi, yn, y_sigma) + dll(dₚ, x_sigma)
        else
            -Inf
        end
        
        # Accept or reject proposal
        if log(rand()) < (llₚ - ll)
            δ ≈ 0 || (jump_sigma = 2.9*abs(δ))
            copyto!(d, dₚ)
            ll = llₚ
        end
        lldist[i] = ll
    end
    for i in 1:nsteps
        # Reset proposals
        copyto!(dₚ, d)

        # Perturb displacement in one region using a Gaussian kernel
        r = rand(1:length(d))
        δ = randn()*jump_sigma
        bandwidth = rand()*x_sigma
        kernel = Normal(xn[r]+d[r], bandwidth)
        for i in eachindex(dₚ)
            dₚ[i] += δ * pdf(kernel, xn[i]+d[i]) * bandwidth
        end
        
        # Apply perturbation
        @. x1ₚ = xn + dₚ

        # Calculate log likelihood
        llₚ = if issorted(x1ₚ)
            linterp1!(ymi, xc, yc, x1ₚ)
            yll(ymi, yn, y_sigma) + dll(dₚ, x_sigma)
        else
            -Inf
        end
        
        # Accept or reject proposal
        if log(rand()) < (llₚ - ll)
            acceptancedist[i] = true
            δ ≈ 0 || (jump_sigma = 2.9*abs(δ))
            copyto!(d, dₚ)
            ll = llₚ
        end

        # Record results
        ddist[:,i] = d
        lldist[i+burnin] = ll
        proposaldist[i] = jump_sigma
    end
    return ddist, lldist, proposaldist, acceptancedist
end

## --- # Generate a noisy initial record
PWT = readdlm("Data/NARR_PWT.csv", ',')  # Read CSV as a matrix
PRECIP = readdlm("Data/NARR_PRECIP.csv", ',')  # Read CSV as a matrix
DATES = readdlm("Data/years.csv", ',')  # Read CSV as a matrix
LEVELS = readdlm("Data/NARR_levels.csv", ',')  # Read CSV as a matrix


dta, header = readdlm("Data/Denali/Isotopes/DEN13A_WithTime.csv", ',', header=true);
ISO = DataFrame(dta, vec(header))
TS1 = readdlm("Data/Denali/DIC1_Timescale_20230620.csv", ',')[2:end,:];
ISO.TopTime = linterp1(TS1[:,2],TS1[:,1],ISO.TopDepth)
ISO.MidTime = linterp1(TS1[:,2],TS1[:,1],ISO.MidDepth)
ISO.BotTime = linterp1(TS1[:,2],TS1[:,1],ISO.BotDepth)
tI = 1980 .<= ISO.MidTime .< 2011
ISO = ISO[tI,:]
ISO = reverse(ISO)

xc = DATES[:,1] .+ (DATES[:,2] .- 1) ./ 26
yc = PWT[:,findfirst(LEVELS .== 500)[2]]

xn = ISO.MidTime
yn = ISO.Deuterium

yc = (yc .- mean(yc)) ./ std(yc)
yn = (yn .- mean(yn)) ./ std(yn)
h = plot(xc,yc, linewidth = 2, label="Central record", framestyle=:box)
plot!(h, xn,yn, label="noisy record")


## --- Run the MCMC
nsteps = 2*10^6
burnin = 10^6
ddist, lldist, proposaldist, acceptancedist = align(xc, yc, xn, yn; x_sigma=1, y_sigma=0.4, nsteps,burnin)


## ------ PLOTS LOG LIKELIHOOD OVER EACH TIMESTEP ------
hl1 = plot(((-burnin+1):nsteps)/1e3,lldist,
    # xlabel = "Step number (1,000s)",
    ylabel = "Log likelihood",
    framestyle = :box,
    label = "",
    xlims = ((-burnin+1),nsteps) ./ 1e3
)
hl2 = plot(((-burnin+1):nsteps)/1e3,lldist,
    xlabel = "Step number (1,000s)",
    ylabel = "Log likelihood",
    framestyle = :box,
    label = "",
    ylims = (-2500, -750),
    xlims = ((-burnin+1),nsteps) ./ 1e3
)
YT,a = yticks(hl1)[1]
yticks!(hl1,YT,string.(Int.(YT)))
YT,a = yticks(hl2)[1]
yticks!(hl2,YT,string.(Int.(YT)))
vspan!(hl1,[(-burnin+1)/1e3,0]; alpha = 0.2, color="black",label="")
vspan!(hl2,[(-burnin+1)/1e3,0]; alpha = 0.2, color="black",label="")
YL = ylims(hl1)
annotate!(hl1,(-burnin+1)/1e3/2,YL[1]+(YL[2]-YL[1])/2,"burn in",fontweight="bold")
#F = text("").font
#F.rotation = 90
#annotate!(hl1,nsteps/1e3/2,YL[1]+(YL[2]-YL[1])/2,text("MCMC",F),fontweight="bold",)
annotate!(hl1,nsteps/1e3/2,YL[1]+(YL[2]-YL[1])/2,"MCMC",fontweight="bold",)
YL = ylims(hl2)
annotate!(hl2,(-burnin+1)/1e3/2,YL[1]+(YL[2]-YL[1])/2,"burn in",fontweight="bold")
#annotate!(hl2,nsteps/1e3/2,YL[1]+(YL[2]-YL[1])/2,text("MCMC",F),fontweight="bold",)
annotate!(hl2,nsteps/1e3/2,YL[1]+(YL[2]-YL[1])/2,"MCMC",fontweight="bold",)
l = @layout [a ; b]
p = plot(hl1,hl2,layout = l,dpi=400,size=[600 400])
#p = plot(hl2,dpi=500,size=[600 300])
savefig(p,"Figures/lldist.png")

## -------------------- Analysis -----------------------
# Calculates mean pertubation for each sample
mu_d = nanmean(ddist, dims=2)

# Loads Denali ice core chemistry data
dta, header = readdlm("Data/Denali/DIC2_IC_ICPMS.csv", ',', header=true); # chemistry data
chem = DataFrame(dta, vec(header))
TS = readdlm("Data/Denali/DIC2_Timescale_20230620.csv", ',')[2:end,:]; # original age scale
chem.TopT = linterp1(TS[:,2],TS[:,1],chem."Top Depth") # interpolates original age scale to ensure accuracy
chem.BotT = linterp1(TS[:,2],TS[:,1],chem."Bot Depth") # interpolates original age scale to ensure accuracy
 
 # Interpolates chemistry data onto isotope timescale
X = 1980:0.01:2011
newMSA = zeros(length(X),1)
newMg = zeros(length(X),1)
newCa = zeros(length(X),1)
for i = 1:size(chem,1)
    topT = chem.TopT[i]
    botT = chem.BotT[i]
    if botT <= 2011 && topT >= 1980
        j,topI = findmin(abs.(X .- topT))
        j,botI = findmin(abs.(X .- botT))
        newMSA[botI:topI] .= chem."IC MSA (ppb)"[i]
        newMg[botI:topI] .= chem."IC Ca2+ (ppb)"[i]
        newCa[botI:topI] .= chem."IC Mg2+ (ppb)"[i]
    end
end
MSAinterp = zeros(size(ISO,1),1)
Mginterp = zeros(size(ISO,1),1)
Cainterp = zeros(size(ISO,1),1)
for i = 1:size(ISO,1)
    topT = ISO.TopTime[i]
    botT = ISO.BotTime[i]
    j,topI = findmin(abs.(X .- topT))
    j,botI = findmin(abs.(X .- botT))
    MSAinterp[i] = mean(newMSA[botI:topI])
    Mginterp[i] = mean(newMg[botI:topI])
    Cainterp[i] = mean(newCa[botI:topI])
end

# Plots chemistry data aligned and unaligned
l = @layout [a ; b; c; d]
p1=plot(xn,yn, lw=1, label="δD",legend=:outerright)
plot!(xn+mu_d,yn, lw=2, label="δD adjusted")
plot!(xc, yc, lw=1, label="PWT 500 hPa", color = "black")
p3=plot(xn,MSAinterp, lw=1, label="MSA",legend=:outerright)
plot!(xn+mu_d,MSAinterp, lw=2, label="MSA adjusted")
p2=plot(xn,Mginterp, lw=1, label="Mg",legend=:outerright)
plot!(xn+mu_d,Mginterp, lw=2, label="Mg adjusted")
p4=plot(xn,Cainterp, lw=1, label="Ca",legend=:outerright)
plot!(xn+mu_d,Cainterp, lw=2, label="Ca adjusted")
XL = (2000,2010)
xlims!(p1,XL)
xlims!(p3,XL)
xlims!(p2,XL)
xlims!(p4,XL)
xticks!(p1,1980:2010)
xticks!(p2,1980:2010)
xticks!(p3,1980:2010)
xticks!(p4,1980:2010)
ylims!(p3,(0,15))
ylims!(p2,(0,40))
ylims!(p4,(0,6))
ylabel!(p1,"norm.")
ylabel!(p3,"ppb")
ylabel!(p2,"ppb")
ylabel!(p4,"ppb")
plot(p1, p2, p3, p4, layout = l, size = [800 550], dpi=500)
savefig("Figures/MC_Method1.png")

#= 
# Plots layer thickness against precipitation
# NOT CORRECTED FOR LAYER THINNING !!!! IGNORE FOR NOW !
xNew = xn + mu_d
precipInterp = zeros(length(xNew)-1,1)
for i = 1:length(xNew)-1
    t = xNew[i+1]
    b = xNew[i]
    j,tI = findmin(abs.(PRECIP[:,1] .- t))
    j,bI = findmin(abs.(PRECIP[:,1] .- b))
    precipInterp[i] = sum(PRECIP[bI:tI,2])
end
precipInterp = precipInterp[:,1]

ISO.Length = ISO.BotDepth .- ISO.TopDepth
acc = ISO.Length[2:end] ./ (12*(xNew[2:end] .- xNew[1:end-1]))
plot(xNew[2:end],acc, lw=2, label="Adjusted chronology",legend=:topleft, size = [700 400], dpi=500)
ylabel!("m month⁻¹")
ylims!((0,2))
p=twinx();
plot!(p,xNew[2:end],precipInterp,lw=0.5,color="black", label="Reanalysis",legend=:topright)
xlims!((2000,2010))
ylims!(p,(0,400))
ylabel!(p,"precip")
xticks!(1980:2010)
savefig("Figures/MC_Method2.png")
cor(acc,precipInterp)
=#


X = ISO.MidDepth;
p1 = plot(dpi=500, yflip = true, xflip = true)
p2 = plot(dpi=500, yflip = true, xflip = true)
for i = eachindex(X)
    d = xn[i] .+ ddist[i,:]
    xMin = floor(minimum(d)*100)/100
    xMax = ceil(maximum(d)*100)/100
    Nbins = 60
    binWidth = (xMax - xMin)/(Nbins - 1)
    xEdges = xMin:binWidth:xMax
    xCenter = (xEdges[1:end-1] .+ xEdges[2:end]) / 2
    Y = histcounts(d,xEdges)
    Y = (Y / length(d)) / binWidth
    plot!(p1,xCenter,X[i] .- Y/100,label="",fillcolor="#4164d9",fillrange = X[i],linecolor="black")
    plot!(p2,xCenter,X[i] .- Y/100,label="",fillcolor="#4164d9",fillrange = X[i],linecolor="black")
end
plot!(p1,ISO.MidTime,ISO.MidDepth,color="red",label="")
ylabel!(p1,"Depth (m)")
xlabel!(p1,"Year")
annotate!(1985,16,"a",fontweight="bold")

plot!(p2,ISO.MidTime,ISO.MidDepth,color="red",label="")
ylims!(p2,(26,36))
XL = (2002,2005)
YL = ylims(p2);
xlims!(p2,XL)
xticks!(p2,1980:2010)
ylabel!(p2,"Depth (m)")
xlabel!(p2,"Year")
annotate!(XL[1]+0.5,YL[1]+0.5,"b",fontweight="bold")
vspan!(p1,[XL[1], XL[2]]; alpha = 0.2, color="black",label="")
l = @layout [a ; b]
p = plot(p1,p2,layout=l,size=[700 500],dpi=500)
savefig(p,"Figures/MC_Agescale.png")



## --- End of File
