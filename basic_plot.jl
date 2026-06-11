include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using Printf

n = 160
a = 2.46

lat_angle = pi/3.0
a1 = [a, 0.0]
a2 = [a*cos(lat_angle), a*sin(lat_angle)]

latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)

# HexUtils.create_honeycomb_lattice!(latA1, latB1, a, a1, a2, false)
HexUtils.create_honeycomb_lattice_2!(latA1, latB1, a1, a2, false)
HexUtils.rotate_lattice!(latA1, pi/6+pi, [0.0, 0.0])
HexUtils.rotate_lattice!(latB1, pi/6+pi, [0.0, 0.0])

a1 = [a*sqrt(3)/2, a/2]
a2 = [a*sqrt(3)/2, -a/2]
d1 = (1/3)*(a1+a2)


ax1 = subplot(111,aspect=1)
ax1.scatter(latA1[:,1], latA1[:,2], s=60, color="blue")
ax1.scatter(latB1[:,1], latB1[:,2], s=60, color="orange")

ax1.quiver(0.0, 0.0, a1[1], a1[2], angles="xy", scale_units="xy", scale=1)
ax1.quiver(0.0, 0.0, a2[1], a2[2], angles="xy", scale_units="xy", scale=1)
ax1.quiver(0.0, 0.0, d1[1], d1[2], angles="xy", scale_units="xy", scale=1)
ax1.quiver(d1[1], 0.0, d1[1], d1[2], angles="xy", scale_units="xy", scale=1)
ax1.quiver(2*d1[1], 0.0, d1[1], d1[2], angles="xy", scale_units="xy", scale=1)

ax1.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
ax1.set_xlim([-4.2, 5.4])
ax1.set_ylim([-4.5, 4.5])

tight_layout()

show()
