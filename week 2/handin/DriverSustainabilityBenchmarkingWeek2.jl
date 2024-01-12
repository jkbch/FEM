# DriverSustainabilityBenchmarkingWeek2.jl

# This is a driver script that will be used for assessment of
# week 2 handins in DTU course 02623 Finite Element Method for Partial
# Differential Equations.
#
# Participants of the course supplies the following drivers.
# Test that they work by running this code before submission)
#   - Driver28b.jl
#   - Driver28c.jl
#
# Developed by Allan P. Engsig-Karup.

using Printf
using LinearAlgebra
using SparseArrays
using SymRCM
using AMD
using Plots

include("Driver28.jl")

# TODO PUT IN YOUR GRROUP NO AND STUDENT IDs FOR THE GROUP HERE
groupNo = 18  # Put your group no here.
groupStudentIDs = "s214722, s214695, s204354"  # Put your student id's here.

# PATHS
cd(@__DIR__)  # change to directory
dirThisScript = pwd()  # store path
dirStoreResults = "/home/jakob/Uni/Fem/code/week 2/handin/plots/" #"/Users/apek/02623/Handinresults/"

# PARAMETERS FOR SUUSTAINABILITY CALCULATION
CO2intensity = 0.285  # [kg CO2/kWh], https://communitiesforfuture.org/collaborate/electricity-map/
PowerEstimate = 60  # [kW]

# PARAMETERS FOR THE SELECTED EXERCISES (DO NOT CHANGE)
# Define input parameters
x0 = 0.0
y0 = 0.0
L1 = 1.0
L2 = 1.0
noelms1 = 40
noelms2 = 50
lam1 = 1.0
lam2 = 1.0
f(x, y) = cos(π * x) * cos(π * y)
qt(x, y) = 2 * π^2 * cos(π * x) * cos(π * y)

# EXECUTE CODE
# Let's call the FEM BVP 2D Solver you produced
# time the code using @time

VX, VY, EToV, U = Driver28b(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, qt)

# Call Group 30 solver
fac = 100
total_time = 0.0
for i in 1:fac
    t = @elapsed begin
        Driver28b(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, qt)
    end
    global total_time += t
end
CPUtime1 = total_time / fac
DOF1 = length(U)

x0 = -1.0
y0 = -1.0
L1 = 2.0
L2 = 2.0

VX2, VY2, EToV2, U2 = Driver28c(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, qt)

total_time = 0.0
for i in 1:fac
    t = @elapsed begin
        Driver28c(x0, y0, L1, L2, noelms1, noelms2, lam1, lam2, f, qt)
    end
    global total_time += t
end
CPUtime2 = total_time / fac
DOF2 = length(U2)

CO2eq1 = CPUtime1 / 3600 * PowerEstimate / 1000 * CO2intensity
CO2eq2 = CPUtime2 / 3600 * PowerEstimate / 1000 * CO2intensity

plot1 = plot(VX, VY, U, st=:surface, xlabel="x", ylabel="y", zlabel="u(x,y)", titlefontsize=5,
    title=@sprintf("2.8b. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e",
    groupStudentIDs, CPUtime1, DOF1, noelms1, noelms2, CO2eq1))

plot2 = plot(VX2, VY2, U2, st=:surface, xlabel="x", ylabel="y", zlabel="u(x,y)", titlefontsize=5,
    title=@sprintf("2.8c. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e",
    groupStudentIDs, CPUtime2, DOF2, noelms1, noelms2, CO2eq2))

display(plot(plot1, plot2, layout=(1,2), size=(1307, 576)))

# STORE THE RESULTS
filename = @sprintf("Week2ResultsGroup%d.txt", groupNo)
fid = open(joinpath(dirStoreResults, filename), "w")  # open file identifier
println(fid, groupNo)
println(fid, groupStudentIDs)
println(fid, @sprintf("%.4e %d %.4e", CPUtime1, DOF1, CO2eq1))
println(fid, @sprintf("%.4e %d %.4e", CPUtime2, DOF2, CO2eq2))
close(fid)  # close file identifier
