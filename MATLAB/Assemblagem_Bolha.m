
clear all
clf

% Criacao da particao [0,1]

N = 10;
dim=2*N-1;
p = linspace(0,1,N);

% Nesta expressao podemos fazer a manipulacao de f%
f = @(x) ( sin(x) );

Ahat = sparse(dim,dim);
bhat = zeros(dim,1);

for k=1:(N-1)
    x1 = p(k);
    x2 = p(k+1);
    %function handle porque vamos estar sempre a avaliar
    phi{:,1} = @(x) (x-x2)./(x1-x2);
    phi{:,2} = @(x) (x-x1)./(x2-x1);
    phi{:,3} = @(x) (x-x1).*(x2-x); %funcao bolha%

    dphi{:,1} = @(x)  1/(x1-x2);
    dphi{:,2} = @(x)  1/(x2-x1);
    dphi{:,3} = @(x) -2*x+x2+x1;

    enum = [k k+1 N+k];
    %Regra de quadratura de Gauss
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


figure(1)
hold on;

%Plot da solucao com bolha
for k=1:N-1
    x1 = p(k);
    x2 = p(k+1);
    phi{:,1} = @(x) (x-x2)./(x1-x2);
    phi{:,2} = @(x) (x-x1)./(x2-x1);
    phi{:,3} = @(x) (x-x1).*(x2-x);
    t = linspace(p(k),p(k+1));
    plot(t,phi{1}(t)*u(k)+phi{2}(t)*u(k+1)+phi{3}(t)*u(k+N));
end
hold off;

figure(2)
hold on;

%Plot da solucao sem bolha
for k=1:N-1
    x1 = p(k);
    x2 = p(k+1);
    phi{:,1} = @(x) (x-x2)./(x1-x2);
    phi{:,2} = @(x) (x-x1)./(x2-x1);
    t = linspace(p(k),p(k+1));
    plot(t,phi{1}(t)*u(k)+phi{2}(t)*u(k+1));
end
ylim([0 0.07]);
hold off;