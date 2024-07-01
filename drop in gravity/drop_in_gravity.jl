using Plots
using DelimitedFiles

# Physics settings
gx = 0.0
gy = 10000*(-9.8)
restRho = 1000
RT = 2000
h = 16
hh = h*h
m = 65
nu = 250
dt = 0.001


# Smoothing kernels 
Wp6 = 315/(64*pi*h^9) # density
Wspiky_grad = -45/(pi*h^6) # pressure 
Wnu_lap = 45/(pi*h^6) # viscosity 


# Wall-collision parameters 
deltaOffset = h
collisionDamping = 0.5

# Domain size
zoneWidth = 1000
zoneHeight = 1000

# Initialize arrays
X = []
Y = []
U = []
V = []
Fx = []
Fy = []
Rho = []
P = []


# Init a fluid vessel
for y in deltaOffset:h:(zoneHeight*0.3 + deltaOffset)
    for x in deltaOffset:h:(zoneWidth - deltaOffset)
        jitter = rand()
        push!(X, x + jitter)
        push!(Y, y)
        push!(U, 0)
        push!(V, 0)
        push!(Fx, 0)
        push!(Fy, 0)
        push!(Rho, 0)
        push!(P, 0)
    end
end

# Init a droplet
x0 = 0.5*zoneWidth
y0 = 0.8*zoneHeight
r0 = 0.1*zoneHeight
for y in deltaOffset:h:(zoneHeight - deltaOffset)
    for x in deltaOffset:h:(zoneWidth - deltaOffset)
        if ((x - x0)^2 + (y - y0)^2 <= r0^2)
            jitter = rand()
            push!(X, x + jitter)
            push!(Y, y)
            push!(U, 0)
            push!(V, 0)
            push!(Fx, 0)
            push!(Fy, 0)
            push!(Rho, 0)
            push!(P, 0)
        end
    end
end

# Get total count of particles
lenp = length(X)

# Render bounds
Xb = [0, 0, zoneWidth, zoneWidth, 0]
Yb = [0, zoneHeight, zoneHeight, 0, 0]
plot(Xb, Yb, seriestype = :shape, linecolor=:black, fillalpha=0, aspect_ratio=:equal, xlabel="X [m]", ylabel="Y [m]")

# Render particles
s = h*2
scatter!(X, Y, markersize=3, markercolor=:blue, aspect_ratio=:equal, legend=false)

num = 10000
tmax = 1

# Time loop
t = 0
while true
    # Update density and pressure
    for ip in 1:lenp
        Rho[ip] = 0
        for jp in 1:lenp
            rij_x = X[jp] - X[ip]
            rij_y = Y[jp] - Y[ip]
            r2 = rij_x^2 + rij_y^2
            r = sqrt(r2)
            if r < h
                Rho[ip] += m * Wp6*(hh - r2)^3
            end
        end
        P[ip] = RT * (Rho[ip] - restRho)
    end

    # Update forces
    for ip in 1:lenp
        fpress_x = 0
        fpress_y = 0
        fvisc_x = 0
        fvisc_y = 0
        fsurf_ten_x = 0
        fsurf_ten_y = 0
        for jp in 1:lenp
            if ip != jp
                rij_x = X[jp] - X[ip]
                rij_y = Y[jp] - Y[ip]
                r2 = rij_x^2 + rij_y^2
                r = sqrt(r2)
                if r < h
                    Pave = (P[ip] + P[jp])/2
                    nablaW = Wspiky_grad*(h - r)^2
                    lapW = Wnu_lap*(h - r)
                    mDivRho = m/Rho[jp]
                    fpress_x -= mDivRho * (rij_x/r) * Pave * nablaW
                    fpress_y -= mDivRho * (rij_y/r) * Pave * nablaW
                    fvisc_x += mDivRho * nu * (U[jp] - U[ip]) * lapW
                    fvisc_y += mDivRho * nu * (V[jp] - V[ip]) * lapW
                end
            end
        end
        fgrav_x = gx * Rho[ip]
        fgrav_y = gy * Rho[ip]
        Fx[ip] = fpress_x + fvisc_x + fgrav_x
        Fy[ip] = fpress_y + fvisc_y + fgrav_y
    end

    # Update velocity and positions
    for ip in 1:lenp
        U[ip] += dt*Fx[ip]/Rho[ip]
        V[ip] += dt*Fy[ip]/Rho[ip]
        X[ip] += dt*U[ip]
        Y[ip] += dt*V[ip]

        # Boundary conditions
        if (X[ip] - deltaOffset) < 0
            U[ip] = -U[ip] * collisionDamping
            X[ip] = deltaOffset
        end
        if (X[ip] + deltaOffset) > zoneWidth
            U[ip] = -U[ip] * collisionDamping
            X[ip] = zoneWidth - deltaOffset
        end
        if (Y[ip] - deltaOffset) < 0
            V[ip] = -V[ip] * collisionDamping
            Y[ip] = deltaOffset
        end
        if (Y[ip] + deltaOffset) > zoneHeight
            V[ip] = -V[ip] * collisionDamping
            Y[ip] = zoneHeight - deltaOffset
        end
    end

    data = hcat(X, Y)  # Склеиваем массивы X и Y
    writedlm(joinpath(desktop_path, "coordinates_$num.txt"), data)  # Сохраняем данные в файл

    num += 1
    t += dt

    if t >= tmax
        break
    end
end
