using Logging, NetCDF
include("particle.jl")

struct HydrodynamicField
  u :: Array{Float32, 3}
  v :: Array{Float32, 3}
  lat :: Array{Float32}
  lon :: Array{Float32}
  time :: Array{Float32}
  invalid_v :: Float32
  invalid_u :: Float32
  lat_max :: Float32
  lat_min :: Float32
  lon_max :: Float32
  lon_min :: Float32
end


function read2DNetcdf(filename :: String)

  if !endswith(filename, ".nc")
    @error "Not a NETCDF file"
  end
  u = ncread(filename, "uo")
  v = ncread(filename, "vo")
  depth = 1
  u = u[:,:,depth,:]
  v = v[:,:,depth,:]
  lat = ncread(filename, "latitude")
  lon = ncread(filename, "longitude")
  time = ncread(filename, "time")
  invalid_u = ncgetatt(filename, "uo", "_FillValue")
  invalid_v = ncgetatt(filename, "vo", "_FillValue")
  lat_max, lat_min = maximum(lat), minimum(lat)
  lon_max, lon_min = maximum(lon), minimum(lon)
  return HydrodynamicField(u,
                           v,
                           lat,
                           lon,
                           time,
                           invalid_v,
                           invalid_u,
                           lat_max,
                           lat_min,
                           lon_max,
                           lon_min)
end

                   