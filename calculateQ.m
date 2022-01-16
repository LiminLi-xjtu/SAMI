
function [Q,Q0] = calculateQ(Y,Q,method)

if strcmp(method,'proc')
    Q = procrust2(Y{1},Y{2} );
    Q0 =0;
elseif strcmp(method,'proc01')
    Q0 = procrust2(Y{1},Y{2} );%由Y计算F得到初
    Q = discreteF(Q0);    
elseif strcmp(method,'descent')
    Q = optimal_F(Q,Y{1}, Y{2}); %更新F
    Q0= 0 ;  
elseif strcmp(method,'descent01')
    Q0 = optimal_F(Q,Y{1}, Y{2}); %更新F
    Q = discreteF(Q0);
elseif strcmp(method,'descent12')
    [Q12, G] = Gsp(Q); %关于F的列稀疏项求导，F12是F的列稀疏值，G是导数
    Q = optimal_F_sp(Q,Y{1}, Y{2}, F12, G); %更新F
    Q0 = 0;
elseif strcmp(method,'union')
    Q = UninCom_F(Q, Y{1}, Y{2});
    Q0 = 0;
end

end