clear all;
N_list = [10 20 40 80 160 320];
f=@(x) (-1).*(x<=1/2)+(1).*(x>1/2);
e=@(x) ((1/2)*x^2-(1/4)*x).*(x<=1/2)+(-(1/2)*(x-1)^2-(1/4)*(x-1)).*(x>1/2);
de=@(x) (x-(1/4)).*(x<=1/2)+(-(x)+1-(1/4)).*(x>1/2);
for i = 1:length(N_list)
    x = linspace(0,1, N_list(i) );
    u = solver1D(f,x);
    [L2error(i),Eerror(i),Herror(i)] = erro_fem1D_linear(x,u,e,de);
end

figure(1)
loglog(1./(N_list-1),L2error,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1).^2,'--or');
legend({'Erro na norma L2','Taxa de h^2'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;

figure(2)
loglog(1./(N_list-1),Eerror,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1),'--or');
legend({'Erro na norma E','Taxa de h'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;

figure(3)
loglog(1./(N_list-1),Herror,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1),'--or');
legend({'Erro na norma H^1','Taxa de h'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;