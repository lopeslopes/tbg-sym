include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf


## finding smallest separation and corresponding angle
angles = Float64[]
AA_sep = Float64[]
BA_sep = Float64[]
AB_sep = Float64[]
BB_sep = Float64[]
AA_vec = []
BA_vec = []
AB_vec = []
BB_vec = []
AA_coord = []
BA_coord = []
AB_coord = []
BB_coord = []

open("AAstack_separations", "r") do f
    for line in eachline(f)
        if startswith(line, "0.")
            push!(angles, parse(Float64, line))
        elseif startswith(line, "AA:")
            aux = split(line, ":")
            push!(AA_sep, parse(Float64, aux[2]))
        elseif startswith(line, "BA:")
            aux = split(line, ":")
            push!(BA_sep, parse(Float64, aux[2]))
        elseif startswith(line, "AB:")
            aux = split(line, ":")
            push!(AB_sep, parse(Float64, aux[2]))
        elseif startswith(line, "BB:")
            aux = split(line, ":")
            push!(BB_sep, parse(Float64, aux[2]))
        elseif startswith(line, "AA_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AA_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BA_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BA_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AB_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AB_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BB_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BB_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AA_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AA_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BA_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BA_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AB_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AB_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BB_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BB_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        end
    end
end

min_index = argmin(AB_sep)
println("Angle:       ", angles[min_index])



## selecting path of minimum separation --
# base_data_path = "data"
# entries = readdir(base_data_path; join=true)
# dataset_dirs = sort(filter(isdir, entries))
# n_datasets = length(dataset_dirs)
#
# path = dataset_dirs[min_index]
## ---------------------------------------

## selecting path mannualy ---------------
AB_stacking = false
angle = angles[min_index]
ang_name = @sprintf("%9.7f", angle)
if AB_stacking
    ang_name = "bernal_"*ang_name
else
    ang_name = "AA_"*ang_name
end
path = "data/"*ang_name*"_2M"
## ---------------------------------------


angle_name = path[findfirst("_",path)[1]+1:findlast("_", path)[1]-1]
angle, moire_period, max_radius, a1_top, a2_top, a1_bot, a2_bot = read_properties(path)

## DISTORTION 1 --------------------------
## changing the a1 and a2 vectors based on the separation vector
## and then generating the whole lattice again using modified lat vectors
cob_matrix = Matrix{Float64}(undef, 2, 2)
cob_matrix[1,1] = a1_top[1]
cob_matrix[1,2] = a2_top[1]
cob_matrix[2,1] = a1_top[2]
cob_matrix[2,2] = a2_top[2]
inv_cob = inv(cob_matrix)

coord_lat_basis = inv_cob * AB_coord[min_index][1:2]
sep_vec_decomp = inv_cob * AB_vec[min_index][1:2]
m_coord = trunc(coord_lat_basis[1])
n_coord = trunc(coord_lat_basis[2])
# m_coord = coord_lat_basis[1]
# n_coord = coord_lat_basis[2]
println("m: ", m_coord)
println("n: ", n_coord)

s1 = sep_vec_decomp[1]/m_coord
s2 = sep_vec_decomp[2]/n_coord
sf = [s1, s2]

# alpha_t1 = [1, 0.0]
# alpha_t2 = [0.0, 1]
# alpha_b1 = [1+s1, 0.0]
# alpha_b2 = [0.0, 1+s2]

alpha_t1 = [1-s1/2, 0.0]
alpha_t2 = [0.0, 1-s2/2]
alpha_b1 = [1+s1/2, 0.0]
alpha_b2 = [0.0, 1+s2/2]

alpha_t1_cart = cob_matrix * alpha_t1
alpha_t2_cart = cob_matrix * alpha_t2
alpha_b1_cart = cob_matrix * alpha_b1
alpha_b2_cart = cob_matrix * alpha_b2

len_alpha_t1 = sqrt(alpha_t1_cart[1]^2 + alpha_t1_cart[2]^2)
len_alpha_t2 = sqrt(alpha_t2_cart[1]^2 + alpha_t2_cart[2]^2)
len_alpha_b1 = sqrt(alpha_b1_cart[1]^2 + alpha_b1_cart[2]^2)
len_alpha_b2 = sqrt(alpha_b2_cart[1]^2 + alpha_b2_cart[2]^2)

println("a1  : ", a1_top)
println("a2  : ", a2_top)

println("\u03B1_t1: ", alpha_t1_cart)
println("\u03B1_t1 length: ", len_alpha_t1)
println("\u03B1_t2: ", alpha_t2_cart)
println("\u03B1_t2 length: ", len_alpha_t2)
println("\u03B1_b1: ", alpha_b1_cart)
println("\u03B1_b1 length: ", len_alpha_b1)
println("\u03B1_b2: ", alpha_b2_cart)
println("\u03B1_b2 length: ", len_alpha_b2)

## ---------------------------------------

## DISTORTION 2 --------------------------
## not distorting for test purposes
# alpha_t1_cart = a1_top
# alpha_t2_cart = a2_top
#
# alpha_b1_cart = a1_bot
# alpha_b2_cart = a2_bot
## ---------------------------------------


## MAKING NEW LATTICES BASED ON THE NEW ALPHA LATTICE VECTORS
n = 2000000

latA1_distorted = zeros(n ÷ 2, 2)
latB1_distorted = zeros(n ÷ 2, 2)
latA2_distorted = zeros(n ÷ 2, 2)
latB2_distorted = zeros(n ÷ 2, 2)

## ROTATING ONE OF THE LATTICES ----------
## method 1: rotating the basis vectors
## then generating the lattice
rotate_point!(alpha_b1_cart, angle, [0.0, 0.0])
rotate_point!(alpha_b2_cart, angle, [0.0, 0.0])

HexUtils.create_honeycomb_lattice!(latA1_distorted, latB1_distorted, alpha_t1_cart, alpha_t2_cart, false)
HexUtils.create_honeycomb_lattice!(latA2_distorted, latB2_distorted, alpha_b1_cart, alpha_b2_cart, false)
## ---------------------------------------
max_radius = maximum(latA1_distorted) - 10.0

## DISTORTION METHOD 3 -------------------
## just to test if the program finds the
## minimum separation point
# latA2_distorted = transpose(transpose(latA2_distorted) .+ AB_vec[min_index][1:2])
# latB2_distorted = transpose(transpose(latB2_distorted) .+ AB_vec[min_index][1:2])
## ---------------------------------------

## ROTATING ONE OF THE LATTICES ----------
## method 2: generating the lattice
## then rotating the whole thing
# HexUtils.create_honeycomb_lattice!(latA1_distorted, latB1_distorted, alpha_t1_cart, alpha_t2_cart, false)
# HexUtils.create_honeycomb_lattice!(latA2_distorted, latB2_distorted, alpha_b1_cart, alpha_b2_cart, false)
#
# rotate_lattice!(latA2_distorted, angle, [0.0, 0.0])
# rotate_lattice!(latB2_distorted, angle, [0.0, 0.0])
## ---------------------------------------


cell_area_bot = cell_area(alpha_b1_cart, alpha_b2_cart)
cell_area_top = cell_area(alpha_t1_cart, alpha_t2_cart)
println("Cell area top lattice: ", cell_area_top)
println("Cell area bot lattice: ", cell_area_bot)

# WRITING POINTS OUT OF OVERLAP
try write_lattice(latA1_distorted, path*"/latticeA1_dist.dat", max_radius)
catch e
    println(e)
end

try write_lattice(latB1_distorted, path*"/latticeB1_dist.dat", max_radius)
catch e
    println(e)
end

try write_lattice(latA2_distorted, path*"/latticeA2_dist.dat", max_radius)
catch e
    println(e)
end

try write_lattice(latB2_distorted, path*"/latticeB2_dist.dat", max_radius)
catch e
    println(e)
end

# finding AA, AB, BA, BB points for the angle and new A2, B2 lattices
treeA1 = KDTree(transpose(latA1_distorted))
treeB1 = KDTree(transpose(latB1_distorted))

tol = 1.0e-3
println("Tolerance: ", tol)
name = @sprintf("%6.4f", tol)

AA = []
BA = []
AB = []
BB = []

n = minimum([size(latA1_distorted)[1], size(latB1_distorted)[1], size(latA2_distorted)[1], size(latB2_distorted)[1]])

for i in 1:n
    indAA, distAA = knn(treeA1, latA2_distorted[i,:], 1)
    indBA, distBA = knn(treeB1, latA2_distorted[i,:], 1)
    indAB, distAB = knn(treeA1, latB2_distorted[i,:], 1)
    indBB, distBB = knn(treeB1, latB2_distorted[i,:], 1)
    if distAA[1] < tol
        push!(AA, latA2_distorted[i,:])
    end
    if distBA[1] < tol
        push!(BA, latA2_distorted[i,:])
    end
    if distAB[1] < tol
        push!(AB, latB2_distorted[i,:])
    end
    if distBB[1] < tol
        push!(BB, latB2_distorted[i,:])
    end
end

latAA = transpose(hcat(AA...))
latBA = transpose(hcat(BA...))
latAB = transpose(hcat(AB...))
latBB = transpose(hcat(BB...))



try write_lattice(latAA, path*"/latticeAA_dist.dat", max_radius)
catch e
    # println(e)
    println("AA lattice is empty!")
end

try write_lattice(latBA, path*"/latticeBA_dist.dat", max_radius)
catch e
    # println(e)
    println("BA lattice is empty!")
end

try write_lattice(latAB, path*"/latticeAB_dist.dat", max_radius)
catch e
    # println(e)
    println("AB lattice is empty!")
end

try write_lattice(latBB, path*"/latticeBB_dist.dat", max_radius)
catch e
    # println(e)
    println("BB lattice is empty!")
end
