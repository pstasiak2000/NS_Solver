#!/usr/bin/env julia
# using AbstractFFTs
using WriteVTK
using FortranFiles
using FFTW
using Mmap
import Printf: @sprintf


struct WaveNumbers
    x::Vector{Float64}
    y::Vector{Float64}
    Nx::Int
    Ny::Int
end

function WaveNumbers(Nx::Int64, Ny::Int64, Lx, Ly)
    kx = rfftfreq(Nx, 2pi * Nx / Lx)
    ky = fftfreq(Ny, 2pi * Ny / Ly)
    return WaveNumbers(kx, ky, length(kx), length(ky))
end

function get_vort!(vx, vy, kk::WaveNumbers)
    vxf = rfft(vx, 1:2)
    vyf = rfft(vy, 1:2)
    wf = similar(vyf)
    for j in 1:kk.Ny, i in 1:kk.Nx
        wf[i,j,1] = im * (kk.x[i] * vyf[i,j,1] - kk.y[j] * vxf[i,j,1])
    end
    w = irfft(wf, kk.Ny, 1:2)
    return w
end


# File containing binary data
folder_fields = "./"

# Output file (without extension)
filename_pvd = "Vel"

# steps to convert
steps = 0:100

# Resolution
Nx =32 
Ny =32 
Nz =32 

# Slab to cut
nslab = 65   


x = range(0, 2pi, length=Nx + 1)[1:end - 1]
y = range(0, 2pi, length=Ny + 1)[1:end - 1]
z = range(0, 2pi, length=Nz + 1)[1:end - 1]

kk = WaveNumbers(Nx, Ny, 2pi, 2pi)

# vx2D = Array{Float64}(undef, Nx, Ny, 1)
# vy2D = Array{Float64}(undef, Nx, Ny, 1)
# vz2D = Array{Float64}(undef, Nx, Ny, 1)
w = Array{Float64}(undef, Nx, Ny, 1)

# data = permutedims(data, (2, 1, 3))

pvd = paraview_collection(folder_fields * filename_pvd)

for it in steps
    itstr = @sprintf "%03d" it
    basename = folder_fields * "V3D.$itstr"
    vtkfile = vtk_grid(basename, x, y, z, compress=false)
        # read!(folder_fields * "vx2D.$itstr.dat", vx2D)
        # read!(folder_fields * "vy2D.$itstr.dat", vy2D)
        # read!(folder_fields * "vz2D.$itstr.dat", vz2D)
        
        # Creat a map to the 3D fields
    uvx = open(folder_fields * "vx.$itstr.dat")
    uvy = open(folder_fields * "vy.$itstr.dat")
    uvz = open(folder_fields * "vz.$itstr.dat")
    Vx = Mmap.mmap(uvx, Array{Float64,3}, (Nx, Ny, Nz))
    Vy = Mmap.mmap(uvy, Array{Float64,3}, (Nx, Ny, Nz))
    Vz = Mmap.mmap(uvz, Array{Float64,3}, (Nx, Ny, Nz))
        # extract the 2D slab
    #vx2D = Vx[:,:,nslab]
    #vy2D = Vy[:,:,nslab]
    #vz2D = Vz[:,:,nslab]


        # Read times
#    timefileUnit = open(read, folder_fields*"tt.$itstr.dat")
#    io = IOBuffer(timefileUnit)
#    tt = read(io,Float64)
#    close(timefileUnit)
#    timefileUnit = FortranFile(folder_fields * "tt.$itstr.dat", "r", access="direct", recl=8)
#    tt = read(timefileUnit,rec=1, Float64)
#    close(timefileUnit)
    tt = it#round(tt, digits=3)
        
        # vtkfile["Time"] = tt
    vtkfile["Velocity"] = (Vx, Vy, Vz)
    #w = get_vort!(vx2D, vy2D, kk)
    #vtkfile["Wz"] = w
    outfiles = vtk_save(vtkfile)
    pvd[tt] = vtkfile

end

vtk_save(pvd)
   
    
