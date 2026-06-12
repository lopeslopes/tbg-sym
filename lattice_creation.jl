include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf

# INITIAL DEFINITIONS
n = 500000

a_top = 2.46
a1_top = [a_top, 0.0]
a2_top = [-a_top*cos(pi/3.0), a_top*sin(pi/3.0)]

AB_stacking = true


# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a1_top, a2_top, false)

max_radius = maximum(latA1) - 10.0

try write_lattice(latA1, "data/latticeA1_500k.dat", max_radius)
catch e
    println(e)
end

try write_lattice(latB1, "data/latticeB1_500k.dat", max_radius)
catch e
    println(e)
end

treeA1 = KDTree(transpose(latA1))
treeB1 = KDTree(transpose(latB1))

# for q in [63.0, 62.0, 61.0, 60.0, 59.0, 58.0, 57.0, 56.0]
for q in [61.0, 60.0, 59.0, 58.0]
    p = 1.0
    # q = 56.0
    angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
    angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
    println(angle_i)
    println(angle_f)
    steps = 100

    for j in 1:steps
        angle = angle_i + j*(angle_f - angle_i)/steps
        ang_name = @sprintf("%9.7f", angle)
        if AB_stacking
            ang_name = "bernal_"*ang_name
        else
            ang_name = "AA_"*ang_name
        end
        create_tbg_system(n, angle, ang_name, treeA1, treeB1, AB_stacking, max_radius)
        a_bot = 2.46
        a1_bot = [a_bot, 0.0]
        a2_bot = [-a_bot*cos(pi/3.0), a_bot*sin(pi/3.0)]
        write_properties(p, q, j, steps, max_radius, a_top, a1_top, a2_top, a_bot, a1_bot, a2_bot, "data/"*ang_name*"_500k/properties.dat")
    end
end
