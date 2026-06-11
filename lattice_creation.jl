include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf


# INITIAL DEFINITIONS
n = 200000

a_top = 2.46
a1_top = [a_top, 0.0]
a2_top = [-a_top*cos(pi/3.0), a_top*sin(pi/3.0)]

AB_stacking = true

rot_axis = [0.0, 0.0]

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a1_top, a2_top, false)

max_radius = maximum(latA1) - 10.0

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
    steps = 200

    for j in 1:steps
        latA2 = zeros(n ÷ 2, 2)
        latB2 = zeros(n ÷ 2, 2)

        a_bot = 2.46
        a1_bot = [a_bot, 0.0]
        a2_bot = [-a_bot*cos(pi/3.0), a_bot*sin(pi/3.0)]

        HexUtils.create_honeycomb_lattice!(latA2, latB2, a1_bot, a2_bot, AB_stacking)

        angle = angle_i + j*(angle_f - angle_i)/steps
        ang_name = @sprintf("%9.7f", angle)
        if AB_stacking
            ang_name = "bernal_"*ang_name
        else
            ang_name = "AA_"*ang_name
        end
        println("Angle in radians: ", angle)
        println("Angle in degrees: ", (angle * 180) / pi)

        # ROTATE SECOND LATTICE BY THE ANGLE
        rotate_lattice!(latA2, angle, rot_axis)
        rotate_lattice!(latB2, angle, rot_axis)

        # tol = 5.0e-3
        # println("Tolerance:        ", tol)
        # name = @sprintf("%6.4f", tol)
        #
        # AA = []
        # BA = []
        # AB = []
        # BB = []
        #
        # for i in 1:div(n,2)
        #     indAA, distAA = knn(treeA1, latA2[i,:], 1)
        #     indBA, distBA = knn(treeB1, latA2[i,:], 1)
        #     indAB, distAB = knn(treeA1, latB2[i,:], 1)
        #     indBB, distBB = knn(treeB1, latB2[i,:], 1)
        #     if distAA[1] < tol
        #         push!(AA, latA2[i,:])
        #     end
        #     if distBA[1] < tol
        #         push!(BA, latA2[i,:])
        #     end
        #     if distAB[1] < tol
        #         push!(AB, latB2[i,:])
        #     end
        #     if distBB[1] < tol
        #         push!(BB, latB2[i,:])
        #     end
        # end
        #
        # latAA = transpose(hcat(AA...))
        # latBA = transpose(hcat(BA...))
        # latAB = transpose(hcat(AB...))
        # latBB = transpose(hcat(BB...))

        try mkdir("data/"*ang_name*"_200k")
        catch e
        end

        # try write_lattice(latAA, "data/"*ang_name*"_200k/latticeAA.dat", max_radius)
        # catch e
        #     println("AA lattice is empty!")
        # end
        #
        # try write_lattice(latBA, "data/"*ang_name*"_200k/latticeBA.dat", max_radius)
        # catch e
        #     println("BA lattice is empty!")
        # end
        #
        # try write_lattice(latAB, "data/"*ang_name*"_200k/latticeAB.dat", max_radius)
        # catch e
        #     println("AB lattice is empty!")
        # end
        #
        # try write_lattice(latBB, "data/"*ang_name*"_200k/latticeBB.dat", max_radius)
        # catch e
        #     println("BB lattice is empty!")
        # end

        # WRITING POINTS OUT OF OVERLAP
        try write_lattice(latA1, "data/"*ang_name*"_200k/latticeA1.dat", max_radius)
        catch e
            println("A1 lat is empty!")
        end

        try write_lattice(latB1, "data/"*ang_name*"_200k/latticeB1.dat", max_radius)
        catch e
            println("B1 lattice is empty!")
        end

        try write_lattice(latA2, "data/"*ang_name*"_200k/latticeA2.dat", max_radius)
        catch e
            println("A2 lat is empty!")
        end

        try write_lattice(latB2, "data/"*ang_name*"_200k/latticeB2.dat", max_radius)
        catch e
            println("B2 lat is empty!")
        end

        write_properties(p, q, j, steps, max_radius, a_top, a1_top, a2_top, a_bot, a1_bot, a2_bot, "data/"*ang_name*"_200k/properties.dat")
    end
end
