using LinearAlgebra, Logging, Random,Distributions


mutable struct Particle
  id :: Int
  lat :: Float32
  lon :: Float32
  grid_x :: Int
  grid_y :: Int
  active :: Bool
end

mutable struct ParticleGroup
  rank :: Int
  loaded_u :: Array{Float32, 3}
  loaded_v :: Array{Float32, 3}
  particle_set :: Vector{Particle}
  group_size :: Int 
  function ParticleGroup(lat₀,
                         lon₀,
                         field)
    particle_set = [Particle(id, lat, lon, 0, 0, true) for (id,(lat, lon)) in enumerate(zip(lat₀, lon₀))]
    rank = 0
    loaded_u = field.u
    loaded_v = field.v
    group_size = length(particle_set)
    new(rank, 
        loaded_u, 
        loaded_v, 
        particle_set,
        group_size)
  end
end

function particle_line(p₁ :: Array{Float64},
                       p₂ :: Array{Float64},
                       n :: Int) 
  m = p₂ - p₁
  ϵ = m/n
  @info "n is $n"
  lat = [p₁[1] + (ϵ[1] * i) for i in 1:n]
  lon = [p₁[2] + (ϵ[2] * i) for i in 1:n]
  return lat, lon
end


function get_normally_distributed_points(a, b, N)
    dist = Normal((a + b) / 2, (b - a) / (2 * sqrt(2)))
    points = rand(dist, N)
    
    return points
end

function random_points(xrange::Vector{Float32}, 
                       yrange::Vector{Float32}, 
                       n::Int)
    xmin, xmax = xrange
    ymin, ymax = yrange
    return [(rand(xmin:xmax), rand(ymin:ymax)) for _ in 1:n]
end
