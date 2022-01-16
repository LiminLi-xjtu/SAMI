clear;
clc

load('data.mat')
X{1} = data_1;
X{2} = data_2;

option.centering =  0;
option.normalize =  0;
option.lapopt = 0;
option.initialP = 'identity';
option.initialY = 'dmap';
option.Qsolver = 'descent';
option.maxit =  1000;
option.epsY = 1e-3;
option.epsobj = 1e-5;


option.dim = 2;
option.alpha = 1;  % alignment
option.beta = 0.001;   % strcture


[Y,Q,P,YY,QQ,PP,obj,cr] = SAMI(X,option);


for i =length(YY):-1:1
   if ~mod(i,10)
    Y = YY{i};
   
    figure,subplot(1,2,1),plot(Y{1}(1,:),Y{1}(2,:),'r*',Y{2}(1,:),Y{2}(2,:),'bo'); hold on
    subplot(1,2,2),scatter(Y{1}(1,:),Y{1}(2,:),40,type1,'filled'); hold on;
    scatter(Y{2}(1,:),Y{2}(2,:),40,type2,'filled'); 
   end
end


