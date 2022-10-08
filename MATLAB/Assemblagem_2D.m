clear all
% funcao source
% Os valores de N indicam o número de vezes a refinar a malha
f = @(x,y)(sin(pi.*x).*sin(pi.*y));

% Criar triangulacao em mesh retangular
N=6;
p = [ 0 10 10 0 ; 0 0 10 10];
t = [ 1 2 3; 1 4 3]';
mesh1=make_mesh(N,p,t);

solver2D(f,N,mesh1);

% Criar triangulacao em mesh exemplo
N=5;
p =[ 0 1 2 0 1 2 ; 0 0 0 1 1 1];
t = [ 1 2 2 3 ; 2 4 3 5 ; 4 5 5 6 ];
mesh2=make_mesh(N,p,t);

solver2D(f,N,mesh2);

% Criar triangulacao em mesh L
N=1;
model = createpde(1);
geometryFromEdges(model,@lshapeg);
lshape=generateMesh(model);
mesh3=make_mesh(N,lshape.Nodes,lshape.Elements);

solver2D(f,N,mesh3);

% Criar triangulacao em mesh Test
N=1;
model = createpde(1);
g = geometryFromEdges(model,@cardg);
testshape=generateMesh(model);
mesh4=make_mesh(N,testshape.Nodes,testshape.Elements);

solver2D(f,N,mesh4); 


