 function [L2error, Eerror] = erro_fem2D(mesh,f,uf,u)
% Devolve as respetivas normas
    % Valores para as normas
    val_L2=0;
    val_E=0;
    % Inicializacao das matrizes 
    P=mesh.p;
    T=mesh.t;
    Nt=size(T,2);

    % Bases referencia 
    phi{:,1} = @(x,y) 1-x-y;
    phi{:,2} = @(x,y) x;
    phi{:,3} = @(x,y) y;
    
    

    for i=1:Nt
        l=0;
        % Funcao erro e derivada do erro
        err = @(x,y) (uf(T(1,i))*phi{1}(x,y)+uf(T(2,i))*phi{2}(x,y)+uf(T(3,i))*phi{3}(x,y))-(u(T(1,i))*phi{1}(x,y)+u(T(2,i))*phi{2}(x,y)+u(T(3,i))*phi{3}(x,y));
        derr = @(x,y) f(x,y)*((uf(T(1,i))*phi{1}(x,y)+uf(T(2,i))*phi{2}(x,y)+uf(T(3,i))*phi{3}(x,y))-(u(T(1,i))*phi{1}(x,y)+u(T(2,i))*phi{2}(x,y)+u(T(3,i))*phi{3}(x,y)));
        % Pesos e pontos
        [t,w] = inttri(20);

        % Matriz At, vetor bt e det 
        At(1,1)=(P(1,T(2,i))-P(1,T(1,i)));
        At(1,2)=(P(1,T(3,i))-P(1,T(1,i)));
        At(2,1)=(P(2,T(2,i))-P(2,T(1,i)));
        At(2,2)=(P(2,T(3,i))-P(2,T(1,i)));
    
        bt(1,1)=P(1,T(1,i));
        bt(2,1)=P(2,T(1,i));

        detAt=det(At);
        
        % Loop de regra de quadratura
        for j=1:length(t)
            Xhat(1,1)=t(1,j);
            Xhat(2,1)=t(2,j);
            Ft=(At*Xhat)+bt;
            val_L2=val_L2+(err(Ft(1),Ft(2)))^2*abs(detAt)*w(j);
            val_E=val_E+(abs(derr(Ft(1),Ft(2))))*abs(detAt)*w(j);
        end
    end
    % Avaliacao das normas
    L2error = sqrt(val_L2);
    Eerror = sqrt(val_E);
end