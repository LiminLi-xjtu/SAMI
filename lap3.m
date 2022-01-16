function L = lap3(Y)

Z=pdist(Y');
S=squareform(Z);

n = size(S,1);
S(find(abs(S)<1e-10)) = 0;
Dis = S;
sigma = median(reshape(Dis,n*n,1));
L = spectral(S,sigma);
end

function L = spectral(W,sigma)
% 谱聚类算法
% 使用Normalized相似变换
% 输入  : W             : N-by-N 矩阵, 即连接矩阵
%        sigma          : 高斯核函数,sigma值
%        num_clusters   : 分类数
%
% 输出  : C : N-by-1矩阵 聚类结果，标签值
%
format long
m = size(W, 1);

%计算相似度矩阵S  相似度矩阵由权值矩阵得到，实践中一般用高斯核函数
W = W.*W;   %平方
W = -W/(2*sigma*sigma);
S = full(spfun(@exp, W));

%获得度矩阵D,此处D为相似度矩阵S中一列元素加起来放到对角线上，得到度矩阵D
D = full(sparse(1:m, 1:m, sum(S)));

% 获得拉普拉斯矩阵 , L =I - D^(-1/2) * S * D^(-1/2)
L = eye(m)-(D^(-1/2) * S * D^(-1/2));

end