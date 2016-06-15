function writeuv(fname, u, v)
    n = length(u)
    f = IOBuffer()#20n)
    for j in 1:n-1
         write(f, string(u[j]))
         write(f, ' ')
         write(f, string(v[j]))
         write(f, ' ')
    end
    write(f, string(u[n]))
    write(f, ' ')
    write(f, string(v[n]))

    open(fname, "w") do g
        write(g, takebuf_string(f))
    end
end
