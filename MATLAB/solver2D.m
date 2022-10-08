function [u] = solver2D(f,N,mesh)
    T=mesh.t;
    P=mesh.p;
    Nt=size(T,2);
    
    bhat=sparse(size(P,2),1);
    Ahat=sparse(size(P,2),size(P,2));
    
    
    for i=1:Nt
        
        % X - Pontos a serem avaliados
        % W - Pesos
        [X,W] = inttri(N);
    
        % Matriz At, vetor bt
    
        At(1,1)=(P(1,T(2,i))-P(1,T(1,i)));
        At(1,2)=(P(1,T(3,i))-P(1,T(1,i)));
        At(2,1)=(P(2,T(2,i))-P(2,T(1,i)));
        At(2,2)=(P(2,T(3,i))-P(2,T(1,i)));
    
        bt(1,1)=P(1,T(1,i));
        bt(2,1)=P(2,T(1,i));
    
        % Calculo det
        detAt=det(At);
    
        % Derivadas
        Nip=size(X,2);
        dL{1} = [ -ones(1,1); -ones(1,1) ];
        dL{2} = [  ones(1,1); zeros(1,1) ];
        dL{3} = [ zeros(1,1);  ones(1,1) ];
    
        for l=1:3
            for j=1:Nip
                % Avaliacao das bases
                phi(:,1) = 1-X(1,j)-X(2,j);
                phi(:,2) = X(1,j);
                phi(:,3) = X(2,j);
                
                % Avaliacao de Ft 
                Xhat(1,1)=X(1,j);
                Xhat(2,1)=X(2,j);
        
                Ft=(At*Xhat)+bt;
    
                bhat(T(l,i))=bhat(T(l,i))+f(Ft(1),Ft(2))*phi(l)*abs(detAt)*W(j);
            end
            for k=1:3
                Ahat(T(l,i),T(k,i))=Ahat(T(l,i),T(k,i))+0.5*transpose(dL{l})*inv(At)*transpose(inv(At))*dL{k}*abs(detAt);
            end
       end
    end

    
    % Pontos de interior e fronteira
    be = find( mesh.e2t(2,:) == 0);
    bind = mesh.edges(:,be);
    bind = unique(bind(:));
    iind = setdiff(1:size(mesh.p,2),bind);

    u =zeros( size(mesh.p,2),1);
    u(iind) = Ahat(iind,iind)\bhat(iind);
    X = mesh.p(1,:);
    Y = mesh.p(2,:);
    t = mesh.t;
    
    figure
    plot_2Dtri_mesh(mesh)
    
    figure
    patch(X(t),Y(t),u(t),u(t));
    view(3);
end