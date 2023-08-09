function ED_S=EffectiveDistance(P,lambda)
%
%% Function EffectiveDistance
%   To compute similarity matrix by using the the effective distance. 
%   Such effective distance is computed based on the sparse reconstruction coefficients.

%% Input parameters:
%  P-        Matrix of data;
%  lambda -  Bandwidth of the expotional function

%% Output parameters:
%  ED_S-     Similarity matrix based on the effective distance 
%
%%  Written by Mingliang Wang  2023.08.09

  [n,size_n]=size(P);

%% To construct the directed graph using sparse representation 

  ED_S = ones(n, n)/n;
  if n>size_n
      option.mode = 'standard';
  else
      option.mode='extend';
  end  
  [sparse_W] = constructW1(P', option);

  mm=sum(sparse_W);
  staticP=zeros(n,n);
  for i=1:size(sparse_W,1)
      staticP(i,:)=sparse_W(i,:)./mm;
  end

  effeDist=1-log(staticP);

  bandWidth=mean(mean(dist(P,P'))); 
  bandRange=lambda*bandWidth;
  ED_S=exp(-(effeDist)/bandRange); 

end

