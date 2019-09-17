function [Av, UpSD, LwSD]=kthNNDcalcRandomNth180830Ripley(FileName2,data10, MaxNeighbor, edges, NoTrial)
%% evaluate by random dist
pnPoints = height(data10);
boxSimu=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];



parfor ii=1:NoTrial
    if ii<NoTrial
        AvAll(:,:,ii)=(NNDfromRandom3DNthPrintRipley(FileName2,pnPoints,boxSimu,edges, MaxNeighbor, 0))';
    elseif ii==NoTrial
        AvAll(:,:,ii)=(NNDfromRandom3DNthPrintRipley(FileName2,pnPoints,boxSimu,edges, MaxNeighbor, 1))';
    end
    
end

Av=mean(AvAll,3);
UpSD=Av+std(AvAll,0,3);
LwSD=Av-std(AvAll,0,3);
Av=Av';
UpSD=UpSD';
LwSD=LwSD';
end