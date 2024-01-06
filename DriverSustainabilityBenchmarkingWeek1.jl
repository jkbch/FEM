using Plots
using LinearAlgebra
using SparseArrays

include("DriverAMR17.jl")

# Identify people behind the work
studentID = "s214722, s214695, s204354"

# Define input parameters (DO NOT CHANGE THIS PART)
funu(x) = exp(-800*(x-0.4)^2) + 0.25*exp(-40*(x-0.8)^2)
func(x) = -1601*exp(-800*(x-0.4)^2) + (-1600*x+640.0)^2*exp(-800*(x-0.4)^2) - 20.25*exp(-40*(x-0.8)^2) + 0.25*(-80*x+64.0)^2*exp(-40*(x-0.8)^2)
x = [0.0, 0.5, 1.0]  # Initial mesh configuration - do not change
M = length(x)
L = 1
c = funu(x[1])
d = funu(x[end])
tol = 1e-4
maxit = 50

# Let's call the FEM BVP 1D Solver with AMR
# time the code using @time

xAMR = []
u = []
iter = 0

# Call Group <X> solver
fac = 100  # we do multiple runs to get the average time
total_time = 0.0
for i in 1:fac
    t = @elapsed begin
        global xAMR, u, iter = DriverAMR17(L, c, d, x, func, tol, maxit)
    end
    global total_time += t
end
CPUtime = total_time / fac
display(CPUtime)

# Plot
p = plot(xAMR, funu.(xAMR), linewidth=2, label="Exact", titlefontsize=7)
scatter!(xAMR, u, markersize=2, label="AMR", xlabel="x", ylabel="u", legend=:topleft)

DOF = length(xAMR)
CO2eq = CPUtime / 3600 * 86 / 1000 * 0.135  # valid for macbook pro (assumed power consumption 105)

title!(string("Group: $studentID, Iter: $iter, Time: ", round(CPUtime, digits=4), " s, DOF: $DOF, CO2e=", round(CO2eq, digits=4), " kg CO2"))
display(p)

if false
    # Element size distribution
    h = diff(xAMR)
    display(histogram(h, xlabel="h", ylabel="# of elements"))
end
