function   [Y,Q,P,YY,QQ,PP,obj,cr] = SAMI(X,option)

%input:
%X size: p time n
%option: 
%option.centering, 1, 0
%option.normalize, 1, 0
%option.initialP, 'pca','identity'
%option.initialY, 'pca','dmap'
%option.alpha
%option.beta: 
%option.Qsovler, 'proc','descent','proc01','descent01'
%option.maxit : number
%option.dim: 
%option.lapopt: 0:'beforepro',1:'afterpro'
%option.eps: for stopping

%output:
% Y: Y{1}, Y{2}, optimal Y with d times n
% Q: optimal Q with n_1 times n_2
% P: optimal P with n_1 times d, n_2 times d

%%% preprocessing
V = length(X);
dim = option.dim;

for v =1:V
   [p(v) n(v)] = size(X{v});
end

if option.lapopt == 0 
   for v = 1:V
      L{v} = lap3(X{v});
   end
end

if option.normalize == 1
  for v = 1: V
    X{v} = normalize0(X{v});
  end
elseif option.normalize == 2 
  for v =1: V  
    X{v} = X{v} - min(X{v}(:));
    X{v} = X{v} / max(X{v}(:));
  end
end


if option.centering == 1
  for v =1: V
    X{v} = centering(X{v});
  end
end



if option.lapopt == 1 
   for v = 1:V
      L{v} = lap3(X{v});
   end
end

 %%% initialization

if strcmp(option.initialP,'pca') && strcmp(option.initialY,'pca')
   for v = 1:V
      [P{v},Y{v}] = pca(X{v}','NumComponents',dim);
   end

end

if strcmp(option.initialY,'dmap')
    for v = 1: V
        Y{v} = diffusionmap(X{v}, dim)';
    end

% 

end

if strcmp(option.initialP,'identity')
    
    for v = 1:V
        P{v} = eye(p(v),dim);
    end
end


[Q,Q0] = calculateQ(Y,zeros(n(1),n(2)),option.Qsolver);


Y_old = Y;
P_old = P;
Q_old = Q;

i = 1;
YY{1} = Y;
PP{1} = P;
QQ{1} = Q;
    
i = i + 1;    
% alternating iterations

alpha = option.alpha;
beta = option.beta;

T = inv((alpha+1)*eye(n(2),n(2))+ beta* L{2} );

while i < option.maxit
    
    %UPDATE P USING Y_old
    for v=1:V
        P{v} = optimal(P_old{v}, Y_old{v}, X{v} ); %algorithm 2 in (Wen, Zaiwe,2013)
    end
    
    
    %UPDATE Y USING Y_old, F_old, P
    
    Y{1} = (P{1}'*X{1} + alpha*Y_old{2}*Q') * inv(alpha* Q_old*Q_old' + beta* L{1} + eye(n(1),n(1)));
    Y{2} = (P{2}'*X{2} + alpha* Y{1}*Q_old ) * T;
    
    
   %update Q
    Q = calculateQ(Y,Q_old,option.Qsolver);
     
    
    %calculate objective value
    
    obj(i,:)=object(X,Y,P,Q,L,alpha,beta,V);
    cr.Q(i) = norm(Q-Q_old,'Fro');
    cr.Y1(i) = norm(Y{1}-Y_old{1},'Fro');
    cr.Y2(i) = norm(Y{2}-Y_old{2},'Fro');
    cr.P1(i) = subspace(P{1},P_old{1});
    cr.P2(i) = subspace(P{2},P_old{2});
    cr.obj(i) = abs(obj(i,4)-obj(i-1,4))/abs(obj(i-1,4));
    
    
   % change new to old
    Y_old = Y;
    P_old = P;
    Q_old = Q;
   
    % save each Y,P,Q
    YY{i} = Y;
    PP{i} = P;
    QQ{i} = Q;
   
   
    % stop 
    if cr.Y1(i) < option.epsY
        break;
    end
    if cr.obj(i) < option.epsobj
        break;
    end
    
    
    i = i + 1
    
end

cr.Y1 = cr.Y1';
cr.Y2 = cr.Y2';
cr.P1 = cr.P1';
cr.P2 = cr.P2';
cr.Q = cr.Q';
cr.obj = cr.obj';

