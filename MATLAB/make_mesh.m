function mesh = make_mesh(N,p,t)
    mesh = inittri(p,t);
    for i=1:N-1
        mesh = refine_tri(mesh);
    end
end