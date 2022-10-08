clear all;
% Script que produz a norma associada ao erro para varias normas e particoes
N_list = [ 3 4 5 6 7];
f = @(x,y)(sin(pi*x)*sin(pi*y));

p = [ 0 1 1 0 ; 0 0 1 1];
t = [ 1 2 3; 1 4 3]';

% Criacao da malha proxima da solucao
meshy=make_mesh(8,p,t);
uf=solver2D(f,7,meshy);

% Loop de comparacao
for i = 1:length(N_list)
    mesh=make_mesh(N_list(i),p,t);
    h(i)=give_h(mesh);
    u=solver2D(f,N_list(i),mesh);
    [L2error(i),Eerror(i)] = erro_fem2D(mesh,f,uf,u);
end

figure(1)
loglog(h,L2error,'k:o');
hold on;
loglog(h,h.^2,'--or');
legend({'Erro na norma L2','Taxa de h^2'},'Location','northwest')
xlabel('Tamanho de h')
hold off;

figure(2)
loglog(h,Eerror,'k:o');
hold on;
loglog(h,h,'--or');
legend({'Erro na norma E','Taxa de h'},'Location','northwest')
xlabel('Tamanho de h')
hold off;