% Codigo adaptado do Professor Antti Hannukainen da Universidade de Aalto [1]

clear all
clf

% Lista das particoes
N_list = [5 10 20 40];

f=@(x)(1+0*x);

t=linspace(0,1,100);

figure(1);
% O comando interpl faz a interpolacao do vetor u da solucao aproximada
for i=1:length(N_list)
    x = linspace(0,1,N_list(i));

    u = solver1D(f,x);

    uh_t=interp1(x,u,t);

    u_t=(1/2)*t.*(1-t);
    
    plot(t,u_t-uh_t);hold on;
end
legend('N=5','N=10','N=20','N=40');