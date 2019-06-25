#!/usr/bin/env julia
using ReadGlobal

function main(args)
    nx,ny,nz,_,_,_ = getdimsize()
    tnx = nx+2
    if testinput(args[1]) && testinput(args[2])

        dtype1,_ = checkinput(args[1],nx,ny,nz)
        u1 = Vector{dtype1}(undef,tnx)
        dtype2,_ = checkinput(args[2],nx,ny,nz)
        u2 = Vector{dtype2}(undef,tnx)
        f1 = open(args[1])
        f2 = open(args[2])
        fo = open(args[3],write=true,append=true)

        nb1 = (dtype1 == Float32 ? 4 : 8)*tnx
        nb2 = (dtype2 == Float32 ? 4 : 8)*tnx

        for _ in 1:(ny*nz)
            unsafe_read(f1,pointer(u1),nb1)
            unsafe_read(f2,pointer(u2),nb2)
            u2 .*= u1

            write(fo,u2)
        end

    else
        error("Strange file size")
    end
end

main(ARGS)