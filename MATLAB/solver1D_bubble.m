function [u] = solver1D_bubble(f,p)
    N=length(p);
    dim=2*N-1;

    Ahat = sparse(dim,dim);
    bhat = zeros(dim,1);

    for k=1:(N-1)
        x1 = p(k);
        x2 = p(k+1);
        phi{:,1} = @(y) (x2-y)./(x2-x1);
        phi{:,2} = @(y) (x1-y)./(x1-x2);
        phi{:,3} = @(x) (x-x1).*(x2-x);
    
        dphi{:,1} = @(x)  1/(x1-x2);
        dphi{:,2} = @(x)  1/(x2-x1);
        dphi{:,3} = @(x) -2*x+x2+x1;
    
        enum = [k k+1 N+k];
        [X,W]=gaussint(2,p(k),p(k+1),1);
    
        for i=1:3
            for l=1:2
                bhat(enum(i))=bhat(enum(i))+W(l)*f(X(l))*phi{i}(X(l));
            end
            for j=1:3
                for l=1:2
                    Ahat(enum(i),enum(j))= Ahat(enum(i),enum(j))+W(l)*dphi{i}(X(l))*dphi{j}(X(l));
                end
            end
        end
    end
    
    idof = setdiff(1:(2*N-1),[1 N]);
    A = Ahat(idof,idof);
    b = bhat(idof,1);
    u(idof) = A\b;
end