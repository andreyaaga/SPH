using Plots
using DelimitedFiles

# Physics settings
gx = 0.0
gy = 0.0
restRho = 1000
RT = 2000
h = 16
hh = h*h
m = 65
nu = 250
dt = 0.001
s_ij = 2
s_ij_stenka = 3

# Smoothing kernels 
Wp6 = 315/(64*pi*h^9) # density 
Wspiky_grad = -45/(pi*h^6) # pressure 
Wnu_lap = 45/(pi*h^6) # viscosity 

# Create a folder on the desktop to save data
desktop_path = joinpath(homedir(), "C:/Users/админ/Desktop/diplom/SimulationData")
mkpath(desktop_path)

# Wall-collision parameters 
deltaOffset = h
collisionDamping = 0.5

# Domain size
zoneWidth = 1500
zoneHeight = 1500

# Initialize arrays
X = []
Y = []
U = []
V = []
Fx = []
Fy = []
Rho = []
P = []

a = 600
b = 400
c = (zoneHeight-a)*0.5-150
d = (zoneWidth-b)*0.5-25

# Init a fluid vessel
for y in ((zoneHeight-a)/2):(h):((zoneHeight+a)/2)
    for x in ((zoneWidth-b)/2):(h):((zoneWidth+b)/2)
        jitter = rand()
        push!(X, x + jitter)
        push!(Y, y + jitter)
        push!(U, 0)
        push!(V, 0)
        push!(Fx, 0)
        push!(Fy, 0)
        push!(Rho, 0)
        push!(P, 0)
    end
end

lenp = length(X)

for y in (c):(h*0.2):(zoneHeight-c)
    for x in (d-5):(h*0.2):(d)
    jitter = rand()
    push!(X, x + jitter)
    push!(Y, y + jitter)
    push!(U, 0)
    push!(V, 0)
    push!(Fx, 0)
    push!(Fy, 0)
    push!(Rho, 0)
    push!(P, 0)
    end
end

for y in (c):(h*0.2):(zoneHeight-c)
    for x in (zoneWidth - d):(h*0.2):(zoneWidth - d + 5)
    jitter = rand()
    push!(X, x + jitter)
    push!(Y, y + jitter)
    push!(U, 0)
    push!(V, 0)
    push!(Fx, 0)
    push!(Fy, 0)
    push!(Rho, 0)
    push!(P, 0)
    end
end



# Get total count of particles
lenp_all = length(X)

# Render bounds
Xb = [0, 0, zoneWidth, zoneWidth, 0]
Yb = [0, zoneHeight, zoneHeight, 0, 0]
plot(Xb, Yb, seriestype = :shape, linecolor=:black, fillalpha=0, aspect_ratio=:equal, xlabel="X [m]", ylabel="Y [m]")

# Render particles
s = h*2
scatter!(X, Y, markersize=3, markercolor=:blue, aspect_ratio=:equal, legend=false)

num = 10000
tmax = 0.5

# Time loop
t = 0
while true
    # Update density and pressure
    for ip in 1:lenp_all
        Rho[ip] = 0
        for jp in 1:lenp_all
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
    for ip in 1:lenp_all
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
                if r < 12*h
                    fsurf_ten_x += s_ij*cos((pi*r)/(24*h))*(rij_x/r)
                    fsurf_ten_y += s_ij*cos((pi*r)/(24*h))*(rij_y/r)
                end
            end
        end

        for jp in (lenp+1):lenp_all
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
                if r < 12*h
                    fsurf_ten_x += s_ij_stenka*cos((pi*r)/(24*h))*(rij_x/r)
                    fsurf_ten_y += s_ij_stenka*cos((pi*r)/(24*h))*(rij_y/r)
                end
            end
        end
        fgrav_x = gx * Rho[ip]
        fgrav_y = gy * Rho[ip]
        Fx[ip] = fpress_x + fvisc_x + fgrav_x + fsurf_ten_x
        Fy[ip] = fpress_y + fvisc_y + fgrav_y + fsurf_ten_y
    end

    # Update velocity and positions
    for ip in 1:lenp
        U[ip] += dt*Fx[ip]/Rho[ip]
        V[ip] += dt*Fy[ip]/Rho[ip]
    end

    for ip in 1:lenp_all
        X[ip] += dt*U[ip]
        Y[ip] += dt*V[ip]
    end

    data = hcat(X, Y)  # Склеиваем массивы X и Y
    writedlm(joinpath(desktop_path, "coordinates_$num.txt"), data)  # Сохраняем данные в файл

    # Render particles
    s = h*1

    # Show time
    num += 1
    t += dt
    title!("t = $t [s], $num")
    display(plot())

    if t >= tmax
        break
    end
end
