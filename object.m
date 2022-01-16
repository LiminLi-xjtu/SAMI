function obj=object(X,Y,P,F,L,beta,gamma,nview)
obj=zeros(1,4);
for i=1:nview
    obj(1)=obj(1) + norm( X{i}-P{i}*Y{i}, 'fro' )^2;
end

obj(2) = obj(2) + norm( Y{1}*F-Y{2}, 'fro' )^2;

for i=1:nview
    obj(3)=obj(3) + trace(Y{i}* L{i} * Y{i}');
end 

obj(4)=obj(1) + beta * obj(2) + gamma * obj(3);
end