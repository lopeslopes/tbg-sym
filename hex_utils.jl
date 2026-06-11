module HexUtils

using LinearAlgebra

export create_honeycomb_lattice!, write_lattice, rotate_lattice!, rotate_point!, read_lattice, read_lattice_3d, magic_angle, write_properties, read_properties, read_properties_raw, cell_area


function create_honeycomb_lattice!(latticeA, latticeB, a1, a2, ab_stacking)
    num_columns = 2*div(isqrt(size(latticeA,1)*2),3)

    v1 = deepcopy(a1)
    v2 = deepcopy(a1)
    rotate_point!(v2, pi/3, [0.0, 0.0])
    v3 = deepcopy(a2)

    d1 = (v1 + v2)*(1/3)

    origin_a = [0.0, 0.0]
    origin_b = deepcopy(d1)

    row = 1
    i = 1
    while(i < size(latticeA,1))
        for j=1:num_columns
            if (i > size(latticeA,1))
                break
            end
            latticeA[i,:] = origin_a + Float64(j-1)*v1
            latticeB[i,:] = origin_b + Float64(j-1)*v1
            i = i + 1
        end
        row  = row + 1
        if (row % 2 == 1)
            origin_a = origin_a + v2
            origin_b = origin_b + v2
        else
            origin_a = origin_a + v3
            origin_b = origin_b + v3
        end
    end

    i0 = div((row ÷ 2) * num_columns + div(num_columns,2), 1)
    lat_origin = latticeA[i0,:]

    latticeA[:,1] = latticeA[:,1] .- lat_origin[1]
    latticeA[:,2] = latticeA[:,2] .- lat_origin[2]
    latticeB[:,1] = latticeB[:,1] .- lat_origin[1]
    latticeB[:,2] = latticeB[:,2] .- lat_origin[2]
    
    if (ab_stacking)
        latticeA[:,1] = latticeA[:,1] .- d1[1]
        latticeA[:,2] = latticeA[:,2] .- d1[2]
        latticeB[:,1] = latticeB[:,1] .- d1[1]
        latticeB[:,2] = latticeB[:,2] .- d1[2]
    end
end

function write_lattice(lattice, filename, max_dist=0.0, min_dist=0.0)
    open(filename, "w") do file
        for i = 1:size(lattice, 1)
            dist = sqrt(lattice[i, 1]^2 + lattice[i, 2]^2)
            if ((max_dist == 0.0) || (dist < max_dist))
                if ((min_dist == 0.0) || (dist > min_dist))
                    println(file, lattice[i, 1], ";", lattice[i, 2])
                end
            end
        end
    end
end

function read_lattice(filename, max_dist=0.0, min_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2])]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    if ((min_dist == 0.0) || (dist > min_dist))
                        push!(lat, aux_v)
                    end
                end
            end
        end
    end
    lat = transpose(hcat(lat...))
    return lat
end

function read_lattice_3d(filename, max_dist=0.0, min_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2]), 0.0]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2 + aux_v[3]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    if ((min_dist == 0.0) || (dist > min_dist))
                        push!(lat, aux_v)
                    end
                end
            end
        end
    end
    lat = transpose(hcat(lat...))
    return lat
end

function rotate_lattice!(lattice, angle, pivot)
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    for i = 1:size(lattice, 1)
        aux1 = lattice[i, :] - pivot
        aux1 = rot_matrix * aux1
        lattice[i, :] .= aux1 .+ pivot
    end
end

function rotate_point!(point, angle, pivot)
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    aux1 = point .- pivot
    aux1 = rot_matrix * aux1
    point .= aux1 .+ pivot
end

function magic_angle(p, q)
    angle = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2))) 
end

function write_properties(p, q, i, steps, max_radius, a_top, a1_top, a2_top, a_bot, a1_bot, a2_bot, filename)
    open(filename, "w") do file
        println(file, "p=", p)
        println(file, "q=", q)
        println(file, "i=", i)
        println(file, "steps=", steps)
        println(file, "max_radius=", max_radius)
        println(file, "a_top=", a_top)
        println(file, "a1_top=", a1_top)
        println(file, "a2_top=", a2_top)
        println(file, "a_bot=", a_bot)
        println(file, "a1_bot=", a1_bot)
        println(file, "a2_bot=", a2_bot)
    end
end

function read_properties(path)
    p = 0
    q = 0
    i = 0
    steps = 0
    max_radius = 0.0
    a_top = 0.0
    a_bot = 0.0
    a1_top = [0.0, 0.0]
    a2_top = [0.0, 0.0]
    a1_bot = [0.0, 0.0]
    a2_bot = [0.0, 0.0]
    open(path*"/properties.dat", "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, "=")
                if aux[1] == "p"
                    aux_p = parse(Float64, aux[2])
                    p = trunc(Int, aux_p)
                elseif aux[1] == "q"
                    aux_q = parse(Float64, aux[2])
                    q = trunc(Int, aux_q)
                elseif aux[1] == "i"
                    aux_i = parse(Float64, aux[2])
                    i = trunc(Int, aux_i)
                elseif aux[1] == "steps"
                    aux_steps = parse(Float64, aux[2])
                    steps = trunc(Int, aux_steps)
                elseif aux[1] == "max_radius"
                    aux_radius = parse(Float64, aux[2])
                    max_radius = aux_radius
                elseif aux[1] == "a_top"
                    aux_atop = parse(Float64, aux[2])
                    a_top = aux_atop
                elseif aux[1] == "a1_top"
                    aux = split(strip(aux[2]), ",")
                    aux_a1top = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a1_top = aux_a1top
                elseif aux[1] == "a2_top"
                    aux = split(strip(aux[2]), ",")
                    aux_a2top = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a2_top = aux_a2top
                elseif aux[1] == "a_bot"
                    aux_abot = parse(Float64, aux[2])
                    a_bot = aux_abot
                elseif aux[1] == "a1_bot"
                    aux = split(strip(aux[2]), ",")
                    aux_a1bot = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a1_bot = aux_a1bot
                elseif aux[1] == "a2_bot"
                    aux = split(strip(aux[2]), ",")
                    aux_a2bot = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a2_bot = aux_a2bot
                end
            end
        end
    end
    angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
    angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
    angle = angle_i + (i/steps)*(angle_f  - angle_i)

    a = 2.46
    moire_period = a/(2*sin(angle/2))

    return angle, moire_period, max_radius, a1_top, a2_top, a1_bot, a2_bot
end

function read_properties_raw(path)
    p = 0
    q = 0
    i = 0
    steps = 0
    max_radius = 0.0
    a_top = 0.0
    a_bot = 0.0
    a1_top = [0.0, 0.0]
    a2_top = [0.0, 0.0]
    a1_bot = [0.0, 0.0]
    a2_bot = [0.0, 0.0]
    open(path*"/properties.dat", "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, "=")
                if aux[1] == "p"
                    aux_p = parse(Float64, aux[2])
                    p = trunc(Int, aux_p)
                elseif aux[1] == "q"
                    aux_q = parse(Float64, aux[2])
                    q = trunc(Int, aux_q)
                elseif aux[1] == "i"
                    aux_i = parse(Float64, aux[2])
                    i = trunc(Int, aux_i)
                elseif aux[1] == "steps"
                    aux_steps = parse(Float64, aux[2])
                    steps = trunc(Int, aux_steps)
                elseif aux[1] == "max_radius"
                    aux_radius = parse(Float64, aux[2])
                    max_radius = aux_radius
                elseif aux[1] == "a_top"
                    aux_atop = parse(Float64, aux[2])
                    a_top = aux_atop
                elseif aux[1] == "a1_top"
                    aux = split(strip(aux[2]), ",")
                    aux_a1top = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a1_top = aux_a1top
                elseif aux[1] == "a2_top"
                    aux = split(strip(aux[2]), ",")
                    aux_a2top = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a2_top = aux_a2top
                elseif aux[1] == "a_bot"
                    aux_abot = parse(Float64, aux[2])
                    a_bot = aux_abot
                elseif aux[1] == "a1_bot"
                    aux = split(strip(aux[2]), ",")
                    aux_a1bot = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a1_bot = aux_a1bot
                elseif aux[1] == "a2_bot"
                    aux = split(strip(aux[2]), ",")
                    aux_a2bot = [parse(Float64, lstrip(aux[1], '[')), parse(Float64, rstrip(aux[2], ']'))]
                    a2_bot = aux_a2bot
                end
            end
        end
    end
    angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
    angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
    angle = angle_i + (i/steps)*(angle_f  - angle_i)

    a = 2.46
    moire_period = a/(2*sin(angle/2))

    return p, q, i, steps, moire_period, max_radius
end

function cell_area(a1, a2)
    a1_3d = [a1[1], a1[2], 0.0]
    a2_3d = [a2[1], a2[2], 0.0]
    area_lat_cell = LinearAlgebra.norm(LinearAlgebra.cross(a1_3d, a2_3d))
    return area_lat_cell
end

end
