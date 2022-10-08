clear all;
%Script que produz a norma associada ao erro para varias normas e particoes
N_list = [10 20 40 80 160 320];
f = @(x)(sin(pi*x));
e = @(x)(sin(pi*x)/(pi^2));
de = @(x)(cos(pi*x)/(pi));
for i = 1:length(N_list)
    x = linspace(0,1, N_list(i) );
    u = solver1D_bubble(f,x);
    [L2error(i),Eerror(i),Herror(i)] = erro_fem1D_bolha(x,u,e,de);
end

figure(1)
loglog(1./(N_list-1),L2error,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1).^3,'--or');
legend({'Erro na norma L2','Taxa de h^3'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;

figure(2)
loglog(1./(N_list-1),Eerror,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1).^2,'--or');
legend({'Erro na norma E','Taxa de h^2'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;

figure(3)
loglog(1./(N_list-1),Herror,'k:o');
hold on;
loglog(1./(N_list-1),1./(N_list-1).^2,'--or');
legend({'Erro na norma H^1','Taxa de h^2'},'Location','northwest')
xlabel('Tamanho da Partição')
hold off;
    