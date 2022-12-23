###################
#FUNCTIONS
###################

using Parsers


# Simple Cubic initial configuration
function sc_init_config(N,L)
       a = L/N^(1/3)
       n = 1
       r = fill(0., N, 3)
       for i in 1:floor(N^(1/3))
               for j in 1:floor(N^(1/3))
                   for k in 1:floor(N^(1/3))
                       r[n,:] = [a*i-L/2 a*j-L/2 a*k-L/2]
                       n += 1
                   end
               end
       end
       return r
end

# Periodic Boundary Conditions
function pbc(rij, L)
       for i in 1:3
              if rij[i] > L/2
                     rij[i] -= L
              elseif rij[i] < -L/2
                     rij[i] += L
              end
       end
       return rij
end
   
   
# Calculate LJ Force and Potential
function find_force_LJ(N,r,L)
       F = fill(0., N, 3)
       V = 0
       for i in 1:N
           for j in i+1:N
              rij = r[i,:]-r[j,:]
              rij = pbc(rij,L)
              d = sqrt(rij[1]^2+rij[2]^2+rij[3]^2)
              if d^2<cutoff^2
                     F[i,:] += (48/d^14-24/d^8)*rij
                     F[j,:] -= (48/d^14-24/d^8)*rij
                     V += 4*(1/d^12-1/d^6) - 4*(1/cutoff^12-1/cutoff^6)
              end
           end
       end
       return F, V
end

# Calculate Kinetic Energy
function kinetic(v)
       T = 0
       for i in 1:N
              T += 1/2*(v[i,1]^2+v[i,2]^2+v[i,3]^2)
       end
       return T
end

# Velocity Verlet time-step
function time_step_vVerlet(N,r,v,L)
       F, V = find_force_LJ(N,r,L)
       for i in 1:N
              r[i,:] = r[i,:] + v[i,:]*dt + 0.5*F[i,:]*dt^2
              r[i,:] = pbc(r[i,:],L)
       end
       v += F * 0.5 * dt
       F, V = find_force_LJ(N,r,L)
       v += F * 0.5 * dt
       return r, v, F, V
end

# Euler time-step
function time_step_Euler_pbc(r,v,L)
       F, V = find_force_LJ(N,r,L)
       for i in 1:N
              r[i,:] = r[i,:] + v[i,:]*dt + 0.5*F[i,:]*dt^2
              r[i,:] = pbc(r[i,:],L)
       end
       v += F*dt
       return r, v, F, V 
end
   
# Andersen thermostat
function therm_Andersen(N,nu,sigma,v)
       for n in 1:N
              if rand() < nu
                     v[n,:] = sigma*randn(3)
              end
       end
end

# Calculate Momentum
function momentum(v)
       p = zeros(1,3)
       for i in 1:N
              p .+= reshape(v[i,:], 1, 3)
       end
       return (p[1]^2+p[2]^2+p[3]^2)/N^2
end

# Initialize velocities with a bimodal distribution
function init_bimodal_velocities(N, T)
       # Initialize the velocities array
       v = zeros(N, 3)
       # Calculate the standard deviation
       sigma = sqrt(T)
       # Loop over all particles
       for i in 1:N
           # Generate a random number
           r = rand()
           # Assign the x velocity
           if r < 0.5
               v[i, 1] = sigma + randn()
           else
               v[i, 1] = -sigma + randn()
           end
           # Generate a random number
           r = rand()
           # Assign the y velocity
           if r < 0.5
               v[i, 2] = sigma + randn()
           else
               v[i, 2] = -sigma + randn()
           end
           # Generate a random number
           r = rand()
           # Assign the z velocity
           if r < 0.5
               v[i, 3] = sigma + randn()
           else
               v[i, 3] = -sigma + randn()
           end
       end
       # Return the velocities array
       return v
end

# Calculate the pressure
function calc_pressure(N, L, Temp, rho, r, F)
       # Calculate the virial
       virial = 0
       for i in 1:N
              ri = r[i,:]
              Fi = F[i,:]
              d = sqrt(ri[1]^2 + ri[2]^2 + ri[3]^2)
              if d^2 < cutoff^2
                     virial += (ri[1]*Fi[1] + ri[2]*Fi[2] + ri[3]*Fi[3])
              end
       end
       # Calculate the pressure using the virial equation
       pressure = rho*Temp + (virial)/(3*L^3)
       return pressure
end


###################
#PARAMETERS
###################

N = 125
rho = 0.05
L = (N/rho)^(1/3)
dt = 0.0001
MCS = Int(1e5)
cutoff = L/2

Temp = 1.2
nu = 0.1
sigma = sqrt(Temp)

###################
#CODE
###################

#initial config
r = sc_init_config(N,L)

#initial velocities
Temp = 100.
v = init_bimodal_velocities(N, Temp)
#v = zeros(N, 3)

#Initial velocity distribution
file3 = open("v_dist_ini.dat", "w")
for i in 1:N
       v_i = v[i,1]
       v_j = v[i,2]
       v_k = v[i,3]
       write(file3, "$v_i \n $v_j \n $v_k \n")
end
close(file3)

file0 = open("trajectory.xyz", "w")
file1 = open("thermodynamics.dat", "w")
#write(file1,"#step     U            Kin        Etot        Temp      Pressure \n")

file2 = open("conservation.dat", "w")
write(file2,"#step  Momentum         Energy \n")

for m in 1:MCS

       print("MCS = ", m, "\r")

       if m>800
              global Temp = 1.2
       end

       #global r,v,F,V = time_step_vVerlet(N,r,v,L)
       global r,v,F,V = time_step_Euler_pbc(r,v,L)
       therm_Andersen(N,nu,sqrt(Temp),v)

       #save trajectory
       if m%(MCS/100) == 0
              #write(file0,"125 \n \n")
              for i in 1:N
                     r_x, r_y, r_z = r[i,1], r[i,2], r[i,3]
                     write(file0,"A $r_x $r_y $r_z \n")
              end
              write(file0,"\n \n")
       end     

       if m>1000
              #save thermodynamics
              if m%(10) == 0

                     Kin = kinetic(v)
                     Kin_aprox = round(Kin; digits=8)

                     Pot = V
                     Pot_aprox = round(Pot; digits=8)

                     Ene = Kin + Pot
                     Ene_aprox = round(Ene; digits=8)

                     Temp_inst = 2*Kin/(3*N-3)
                     Temp_inst_aprox = round(Temp_inst; digits=8)

                     P = calc_pressure(N, L, Temp_inst, rho, r, F)
                     P_aprox = round(P; digits=8)

                     write(file1,"$m,$Pot_aprox,$Kin_aprox,$Ene_aprox,$Temp_inst_aprox,$P_aprox \n")
                     
                     mom = momentum(v)
                     mom_aprox = round(mom; digits=8)

                     write(file2,"$m $mom $Ene \n")
              end 
       end
end

close(file0)
close(file1)
close(file2)

#Final velocity distribution
file4 = open("v_dist_fin.dat", "w")
for i in 1:N
       v_i = v[i,1]
       v_j = v[i,2]
       v_k = v[i,3]
       write(file4, "$v_i \n $v_j \n $v_k \n")
end
close(file4)

#run(`gnuplot -e L=$L final_config.gnu`)
run(`gnuplot -e L=$L traject.gnu`)