# Julia Script to plot the waveform

using Plots
using DelimitedFiles
using FFTW
gr()

x = readdlm("waveform.dat", Float64)

waveform = plot(x[:,1],x[:,2])
