function Scatter_UE2BS = genScatterUEBS(K, S, N, BS, UE, num_scatter_UE2BS,...                               
                PathInCluster,...
                SizeOfScatter)
% UE2BS scatterer
if num_scatter_UE2BS ~= 0        
    Scatter_UE2BS = GenScatterXY(N, UE, BS, num_scatter_UE2BS);      
else
    Scatter_UE2BS.x = [];   
    Scatter_UE2BS.y = [];
end
% add cluster structure
if (num_scatter_UE2BS) ~= 0  
    Scatter_UE2BS_tmp.x = Scatter_UE2BS.x;
    Scatter_UE2BS_tmp.y = Scatter_UE2BS.y;
    Scatter_UE2BS.x = [];
    Scatter_UE2BS.y = [];
    % add cluster structure to the generated Scatter_UE2BS
    for i = 1:PathInCluster
        Scatter_UE2BS.x = [Scatter_UE2BS.x, Scatter_UE2BS_tmp.x-SizeOfScatter/2+SizeOfScatter/PathInCluster*i];
        Scatter_UE2BS.y = [Scatter_UE2BS.y, Scatter_UE2BS_tmp.y];
    end
end
Scatter_UE2BS.theta = atan((Scatter_UE2BS.y-BS.y)./(Scatter_UE2BS.x-BS.x));
% row number equals to the BS antenna dimension
% column number equals to the total path number of scatterers (contain the cluster structure in each scatterer)
Scatter_UE2BS.r = distance(repmat(UE.x, 1, length(Scatter_UE2BS.x))...
                , repmat(UE.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y)...
                ...
                + distance(repmat(BS.x, 1, length(Scatter_UE2BS.x))...
                , repmat(BS.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y);
% Scatter_UE2BS.r = distance(repmat(UE.x, N, length(Scatter_UE2BS.x))...
%                 , repmat(UE.y, N, length(Scatter_UE2BS.y))...
%                 , repmat(Scatter_UE2BS.x, N, 1)...
%                 , repmat(Scatter_UE2BS.y, N, 1))...
%                 ...
%                 + distance(repmat(BS.xAll, 1, length(Scatter_UE2BS.x))...
%                 , repmat(BS.yAll, 1, length(Scatter_UE2BS.y))...
%                 , repmat(Scatter_UE2BS.x, N, 1)...
%                 , repmat(Scatter_UE2BS.y, N, 1));
Scatter_UE2BS.r_LastHop = distance(repmat(BS.x, 1, length(Scatter_UE2BS.x))...
                , repmat(BS.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x ...
                ,Scatter_UE2BS.y);
% Scatter_UE2BS.alpha = 1./distance(repmat(UE.x, N, length(Scatter_UE2BS.x))...
%                 , repmat(UE.y, N, length(Scatter_UE2BS.y))...
%                 , repmat(Scatter_UE2BS.x, N, 1)...
%                 , repmat(Scatter_UE2BS.y, N, 1))...
%                 ...
%                 .*1./distance(repmat(BS.x, N, length(Scatter_UE2BS.x))...
%                 , repmat(BS.y, N, length(Scatter_UE2BS.y))...
%                 , repmat(Scatter_UE2BS.x, N, 1)...
%                 , repmat(Scatter_UE2BS.y, N, 1))...
%                 ./sqrt(PathInCluster)...
%                 ...
%                 ./( (4*pi)^(3/2) );
if isempty(Scatter_UE2BS.x) == 1
    Scatter_UE2BS.alpha = [];
else
    Scatter_UE2BS.alpha = 1./distance(repmat(UE.x, N, length(Scatter_UE2BS.x))...
                    , repmat(UE.y, N, length(Scatter_UE2BS.y))...
                    , repmat(Scatter_UE2BS.x, N, 1)...
                    , repmat(Scatter_UE2BS.y, N, 1))...
                    ...
                    .*1./distance(repmat(BS.xAll, 1, length(Scatter_UE2BS.x))...
                    , repmat(BS.yAll, 1, length(Scatter_UE2BS.y))...
                    , repmat(Scatter_UE2BS.x, N, 1)...
                    , repmat(Scatter_UE2BS.y, N, 1))...
                    ./sqrt(PathInCluster)...
                    ...
                    ./(4*pi)^(3/2)...
                    .*S...
                    ...
                    .*(10.^(K.*(distance(repmat(UE.x, N, length(Scatter_UE2BS.x))...
                    , repmat(UE.y, N, length(Scatter_UE2BS.y))...
                    , repmat(Scatter_UE2BS.x, N, 1)...
                    , repmat(Scatter_UE2BS.y, N, 1)) ...
                    + distance(repmat(BS.xAll, 1, length(Scatter_UE2BS.x))...
                    , repmat(BS.yAll, 1, length(Scatter_UE2BS.y))...
                    , repmat(Scatter_UE2BS.x, N, 1)...
                    , repmat(Scatter_UE2BS.y, N, 1))...
                    )/10));

%     Scatter_UE2BS.alpha = 1./distance(repmat(UE.x, N, length(Scatter_UE2BS.x))...
%                     , repmat(UE.y, N, length(Scatter_UE2BS.y))...
%                     , repmat(Scatter_UE2BS.x, N, 1)...
%                     , repmat(Scatter_UE2BS.y, N, 1))...
%                     ...
%                     .*1./distance(repmat(BS.xAll, 1, length(Scatter_UE2BS.x))...
%                     , repmat(BS.yAll, 1, length(Scatter_UE2BS.y))...
%                     , repmat(Scatter_UE2BS.x, N, 1)...
%                     , repmat(Scatter_UE2BS.y, N, 1))...
%                     ./sqrt(PathInCluster)...
%                     ...
%                     ./( (4*pi)^(3/2) );
end

end