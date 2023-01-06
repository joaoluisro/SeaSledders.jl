using Random, Plots, PlotlyJS

function find_interval(v, element)
  v_len = length(v)
  index = argmin([abs(element - v[i]) + abs(element - v[i - 1]) for i in 2:v_len])
  return((index, index + 1))
end

# interpolate in space
function bilinear2DField(x,
                         y,
                         p :: Particle,
                         grid :: Vector{Tuple{Int64, Int64}},
                         fieldData :: Array{Float32, 2},
                         invalid :: Float32)

  grid_x, grid_y = grid[1], grid[2]

  out_of_bounds = (grid_x[2] > size(fieldData)[2] || grid_y[2] > size(fieldData)[1])
  if out_of_bounds
    p.active = false
    return (0.0)
  end 
  v_grid = (fieldData[grid_x[1], grid_y[1]],
            fieldData[grid_x[2], grid_y[1]],
            fieldData[grid_x[1], grid_y[2]],
            fieldData[grid_x[2], grid_y[2]])
  
  if (invalid in v_grid) 
    return (0.0)
  end
  return (v_grid[1] * (fieldData[grid_x[2]] - x) * (fieldData[grid_y[2]] - y)) +
         (v_grid[2] * (x - fieldData[grid_x[1]]) * (fieldData[grid_y[2]] - y)) +
         (v_grid[3] * (fieldData[grid_x[2]] - x) * (y - fieldData[grid_y[1]])) +
         (v_grid[4] * (x - fieldData[grid_x[1]]) * (y - fieldData[grid_y[1]]))
end

function velocity(field,
                  p,
                  x,
                  y,
                  t)

  # interpolate in time
  time = argmin([time - t for time in field.time])
  # interpolate in space
  grid_x, grid_y = find_interval(field.lat,x), find_interval(field.lon, y)
  grid = collect((grid_x, grid_y))
  vx = bilinear2DField(x, y, p, grid, field.u[:, :, time], field.invalid_u)
  vy = bilinear2DField(x, y, p, grid, field.v[:, :, time], field.invalid_v)
  return (vx, vy)
end

function advect(field :: HydrodynamicField,
                particles :: ParticleGroup,
                dₜ :: Float64;
                n_iter :: Int = 100,
                output :: String = "out.csv")
  
  trajectory = []
  for i in 1:n_iter
    t = i * dₜ
    for p in particles.particle_set
      xₜ, yₜ = p.lat, p.lon
      
      uₜ, vₜ = velocity(field, p, xₜ, yₜ, t)
      if(!p.active)
        continue
      end
      xₜ₊₁ = xₜ + uₜ * dₜ 
      yₜ₊₁ = yₜ + vₜ * dₜ
      @info "lat:$xₜ lon: $yₜ uₜ: $uₜ vₜ: $vₜ"

      xₜ = xₜ₊₁
      yₜ = yₜ₊₁

      push!(trajectory, [p.id, i, xₜ, yₜ])
    end
  end 
  @info "trajectory" trajectory
  return trajectory
end

function Animate(field,
                 n_iter,
                 result)

  animation = @animate for i in 1:n_iter
    Plots.scatter(xlims=(minimum(field.lat), maximum(field.lat)),
                  ylims=(minimum(field.lon), maximum(field.lon)),
                  [result[i][1]],
                  [result[i][2]],
                  legend=false)
  end

  return animation
end