include("hex_utils.jl")
using .HexUtils
using PyCall
using PyPlot
using LinearAlgebra
using NearestNeighbors
using Printf
using Combinatorics
using StatsBase
using Crystalline


sg = spacegroup(191, 3)
db = directbasis(191, 3)

base_data_path = "data"
entries = readdir(base_data_path; join=true)
dataset_dirs = sort(filter(isdir, entries))
n_datasets = length(dataset_dirs)

for path in dataset_dirs
    if (path[6:13] != "AA_stack")
        continue
    end
    println(path)
    angle, moire_period, max_radius = read_properties(path)

    steps = 500
    r_step = (max_radius - 0.2)/steps
    for i in 1:500
        radius = 0.2 + i*((max_radius-r_step)-0.2)/steps

        latA1 = transpose(read_lattice_3d(path*"/latticeA1.dat", radius+(r_step/2), radius-(r_step/2)))
        latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat", radius+(r_step/2), radius-(r_step/2)))
        latA2 = transpose(read_lattice_3d(path*"/latticeA2.dat", radius+(r_step/2), radius-(r_step/2)))
        latB2 = transpose(read_lattice_3d(path*"/latticeB2.dat", radius+(r_step/2), radius-(r_step/2)))

        if (isempty(latA1) || isempty(latB1) || isempty(latA2) || isempty(latB2))
            continue
        end
        treeA1 = KDTree(latA1)
        treeB1 = KDTree(latB1)
        treeA2 = KDTree(latA2)
        treeB2 = KDTree(latB2)

        tol = 1e-8
        notable_points = []
        for test_point in eachcol(latA1)
            ops_A1 = []
            ops_B1 = []
            ops_A2 = []
            ops_B2 = []
            ops_no_match = []
            for op in sg
                flag = 0
                gen_rot = cartesianize(op, db)
                new_p = gen_rot.rotation * test_point

                indA1, distA1 = nn(treeA1, new_p)
                indB1, distB1 = nn(treeB1, new_p)
                indA2, distA2 = nn(treeA2, new_p)
                indB2, distB2 = nn(treeB2, new_p)

                if distA1 < tol
                    push!(ops_A1, op)
                    flag = flag + 1
                end
                if distB1 < tol
                    push!(ops_B1, op)
                    flag = flag + 1
                end
                if distA2 < tol
                    push!(ops_A2, op)
                    flag = flag + 1
                end
                if distB2 < tol
                    push!(ops_B2, op)
                    flag = flag + 1
                end
                if flag == 0
                    push!(ops_no_match, op)
                end
            end

            num_op_not_a1_a1 = size(ops_B1)[1] + size(ops_A2)[1] + size(ops_B2)[1] + size(ops_no_match)[1]
            if num_op_not_a1_a1 > 0
                push!(notable_points, test_point)
            end
        end
        notable_points = hcat(notable_points...)

        if isempty(notable_points)
            # println("   Num. notable points: 0")
        else
            println("   Radius: ", radius)
            println("       Num. notable points: ", size(notable_points)[2])
        end
    end
end

# latA1 = transpose(read_lattice_3d(path*"/latticeAB.dat"))
# # latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat", radius+(r_step/2), radius-(r_step/2)))
# # latA2 = transpose(read_lattice_3d(path*"/latticeA2.dat", radius+(r_step/2), radius-(r_step/2)))
# # latB2 = transpose(read_lattice_3d(path*"/latticeB2.dat", radius+(r_step/2), radius-(r_step/2)))
#
# ax1 = subplot(111)
# ax1.scatter(latA1[1, :], latA1[2, :])
#
# ax1.set_aspect("equal")
# show()
