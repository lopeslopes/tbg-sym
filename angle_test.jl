include("hex_utils.jl")
using .HexUtils

function result(l, m, n)
    aux = m^2 + 3*n^2 + 3*m*n + m + 2*n + 1/3 - 3*l^2
    return aux
end 

for l in 1:5000
    for m in 1:5000
        for n in 1:5000
            aux = result(l,m,n)
            if abs(aux) < 0.5
                println(l, ", ", m, ", ", n, ", ", aux)
            end
        end
    end
end

