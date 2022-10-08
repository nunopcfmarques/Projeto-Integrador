% Codigo adaptado do Professor Antti Hannukainen da Universidade de Aalto [1]

clear all
clf
% Criacao da particao no intervalo [0,1]
N = 10;
x = linspace(0,1,N);

% Nesta expressao podemos fazer a manipulacao de f%
f = @(x)(1+0*x);
uf= @(x)(1/2*(1-x)*x);

% Foi usada a funcao sparse porque Ahat tem muitos valores nulos
Ahat = sparse(N,N);
bhat = zeros(N,1);

% Assemblagem

for k=1:(N-1)
    
    % inicio e fim do intervalo
    x1 = x(k);
    x2 = x(k+1);
    
    
    % termos
    phi{:,1} = @(x) (x-x2)./(x1-x2);
    phi{:,2} = @(x) (x-x1)./(x2-x1);


    % valores da derivada de phi. Obtido logo por derivacao de phi em x
    dphi(1) = 1/(x1-x2);
    dphi(2) = 1/(x2-x1);
    
    enum =[k k+1];
    [X,W]=gaussint(2,x(k),x(k+1),1);
    for i=1:2
        for l=1:2
            bhat(enum(i))=bhat(enum(i)) + W(l)*f(X(l))*phi{i}(X(l));
        end
        for j=1:2
            int=0;
            for l=1:2
                int=int + W(l)*dphi(i)*dphi(j);
            end
            Ahat(enum(i),enum(j)) = Ahat(enum(i),enum(j)) + int;
        end
    end
end

    % Sao removidos o primeiro e o ultimo ponto do intervalo
    A=Ahat(2:(N-1),2:(N-1));
    b=bhat(2:(N-1),1);
    
    % Construcao da solucao u
    u(1,1)=0;
    u(N,1)=0;
    u(2:(N-1),1)=A\b;
    
    % Figura da solucao
    figure(1)
    plot(x,u,'r*-')
    figure(2)
    fplot(uf,[0,1]);
    ylim([0 0.14])