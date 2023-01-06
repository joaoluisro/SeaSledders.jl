module SeaSledders

include("field.jl")
include("advection.jl")
include("particle.jl")

filename = "south-atlantic.nc"
field = read2DNetcdf(filename)

#startPoint = [-47.501, -35.361]
#endPoint = [-53.431, -36.527]
#lat, lon = particle_line(startPoint, endPoint, 10)
#@info lat
#@info lon
lat = get_normally_distributed_points(field.lat_min, field.lat_max, 100)
lon = get_normally_distributed_points(field.lon_min, field.lon_max, 100)
pgroup = ParticleGroup(lat, lon, field)

#@info "group" pgroup.particle_set
trajectory = advect(field, pgroup, 0.01)
#lat, lon = [-30.609, -40.989]
#points = random_points([field.lat_min, field.lat_max], [field.lon_min, field.lon_max], 100)
#@info "$points"
#trajectory = advect(field, pgroup, 0.01)
#anim = Animate(field, 100, trajectory)
#gif(anim, "out.gif", fps=30)
end # module SeaSledders