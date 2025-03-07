## Project 4
## EARS 202
## Jacob Chalif
using Pkg
Pkg.activate(".")
Pkg.add(["SignalAlignment","Distances","DelimitedFiles","Statistics","Plots"])
using  Plots,  DelimitedFiles, Distances, SignalAlignment, Statistics

# Reading in data
ISOraw = readdlm("Data/Denali/Isotopes/DEN13A_WithTime.csv", ',')  # Read CSV as a matrix
ISOraw = ISOraw[2:end,:]
time = ISOraw[:,7]
o18 = ISOraw[:,4]
depth = ISOraw[:,3]
ISO = readdlm("Data/ISO_WEEKLY_INTEPROLATED.csv", ',')  # Read CSV as a matrix
PWT = readdlm("Data/NARR_PWT.csv", ',')  # Read CSV as a matrix
T = readdlm("Data/NARR_T.csv", ',')  # Read CSV as a matrix
LEVELS = readdlm("Data/NARR_levels.csv", ',')  # Read CSV as a matrix
DATES = readdlm("Data/years.csv", ',')  # Read CSV as a matrix


## Signal alignment
lvl = 500 # pressure level to use
lvlI = findfirst(LEVELS .== lvl)[2] # select level
pwtNorm = (PWT[:,lvlI] .- mean(PWT[:,lvlI])) ./ std(PWT[:,lvlI]) # normalized PWT
isoNorm = (ISO[:,1] .- mean(ISO[:,1])) ./ std(ISO[:,1]) # normalized water isotopes
signals = [pwtNorm, isoNorm]  # A vector of signals we want to align

master = Index(1)       # Indicate which signal is the master to which the others are aligned
method = Warp(warp_method=DTW(radius=10))
aligned_signals = align_signals(signals, method; master, output=Signals())
aligned_indices = align_signals(signals, method; master, output=Indices())

dateT = DATES[:,1] .+ (DATES[:,2] .- 1) ./ 26  # time vector in decimal years

# plots unaligned and aligned water isotopes
plot(dateT,signals, label=["PWT" "ISO"], l=(:dash, ))
plot!(dateT,[pwtNorm[aligned_indices[1]] isoNorm[aligned_indices[2]]], label=["PWT aligned" "ISO aligned"], c=(1:2)', size=(600, 400))

# saves aligned signals for plotting in Matlab
dta2sav =[dateT pwtNorm isoNorm pwtNorm[aligned_indices[1]] isoNorm[aligned_indices[2]]];
writedlm("AlignedSignals.csv",dta2sav)

# saves aligned indices only, for use with plotting other chemistry in Matlab
writedlm("AlignedIndices.csv",[dateT aligned_indices[2]])
