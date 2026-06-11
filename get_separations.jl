include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf


base_data_path = "data"
entries = readdir(base_data_path; join=true)
dataset_dirs = sort(filter(isdir, entries))
n_datasets = length(dataset_dirs)

smallest_AA_separation = 10000.0
AA_vector_separation = [0.0, 0.0]
AA_point_coord = [0.0, 0.0]

smallest_BA_separation = 10000.0
BA_vector_separation = [0.0, 0.0]
BA_point_coord = [0.0, 0.0]

smallest_AB_separation = 10000.0
AB_vector_separation = [0.0, 0.0]
AB_point_coord = [0.0, 0.0]

smallest_BB_separation = 10000.0
BB_vector_separation = [0.0, 0.0]
BB_point_coord = [0.0, 0.0]

# path = last(dataset_dirs)
for path in dataset_dirs
    angle_name = path[findfirst("/",path)[1]+1:findlast("_", path)[1]-1]
    println(angle_name)

    latA1 = transpose(read_lattice_3d(path*"/latticeA1.dat"))
    latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat"))
    latA2 = transpose(read_lattice_3d(path*"/latticeA2.dat"))
    latB2 = transpose(read_lattice_3d(path*"/latticeB2.dat"))

    treeA1 = KDTree(latA1)
    treeB1 = KDTree(latB1)
    treeA2 = KDTree(latA2)
    treeB2 = KDTree(latB2)

    latA1 = transpose(latA1)
    latB1 = transpose(latB1)
    latA2 = transpose(latA2)
    latB2 = transpose(latB2)

    AA = []
    BA = []
    AB = []
    BB = []

    n = minimum([size(latA1)[1], size(latB1)[1], size(latA2)[1], size(latB2)[1]])

    global smallest_AA_separation = 10000.0
    global AA_vector_separation = [0.0, 0.0]
    global AA_point_coord = [0.0, 0.0]

    global smallest_BA_separation = 10000.0
    global BA_vector_separation = [0.0, 0.0]
    global BA_point_coord = [0.0, 0.0]

    global smallest_AB_separation = 10000.0
    global AB_vector_separation = [0.0, 0.0]
    global AB_point_coord = [0.0, 0.0]

    global smallest_BB_separation = 10000.0
    global BB_vector_separation = [0.0, 0.0]
    global BB_point_coord = [0.0, 0.0]

    for i in 1:div(n,2)
        indAA, distAA = knn(treeA1, latA2[i,:], 1)
        if distAA[1] < smallest_AA_separation
            global smallest_AA_separation = distAA[1]
            global AA_vector_separation = transpose(latA1[indAA,1:3]) - latA2[i,:]
            global AA_vector_separation = AA_vector_separation[1:3]
            global AA_point_coord = latA2[i,:]
        end

        indBA, distBA = knn(treeB1, latA2[i,:], 1)
        if distBA[1] < smallest_BA_separation
            global smallest_BA_separation = distBA[1]
            global BA_vector_separation = transpose(latB1[indBA,1:3]) - latA2[i,:]
            global BA_vector_separation = BA_vector_separation[1:3]
            global BA_point_coord = latA2[i,:]
        end

        indAB, distAB = knn(treeA1, latB2[i,:], 1)
        if distAB[1] < smallest_AB_separation
            global smallest_AB_separation = distAB[1]
            global AB_vector_separation = transpose(latA1[indAB,1:3]) - latB2[i,:]
            global AB_vector_separation = AB_vector_separation[1:3]
            global AB_point_coord = latB2[i,:]
        end

        indBB, distBB = knn(treeB1, latB2[i,:], 1)
        if distBB[1] < smallest_BB_separation
            global smallest_BB_separation = distBB[1]
            global BB_vector_separation = transpose(latB1[indBB,1:3]) - latB2[i,:]
            global BB_vector_separation = BB_vector_separation[1:3]
            global BB_point_coord = latB2[i,:]
        end
    end

    println("AA:        ", smallest_AA_separation)
    println("AA_vec:    ", AA_vector_separation)
    println("AA_coord:  ", AA_point_coord)
    println("BA:        ", smallest_BA_separation)
    println("BA_vec:    ", BA_vector_separation)
    println("BA_coord:  ", BA_point_coord)
    println("AB:        ", smallest_AB_separation)
    println("AB_vec:    ", AB_vector_separation)
    println("AB_coord:  ", AB_point_coord)
    println("BB:        ", smallest_BB_separation)
    println("BB_vec:    ", BB_vector_separation)
    println("BB_coord:  ", BB_point_coord)

    println("")
end
