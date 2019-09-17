function [Av, UpSD, LwSD]=kthNNDcalcRandomDimerBlinkNth180918Ripley(FileName2, data10, MaxNeighbor, edges, NoTrial,ObservationRate, Dist1stN,Dist2ndN,DAngle, PopulationVector, Accu)
% data10: table of the original data
% pnPoints: No of data
% box=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];
% ObserbatioRate: % observation in the lattice (0-1)
% tiltAngleX: degree of rotation of the lattice relatice to the xy plane (0-90)
% tiltAngleZ: degree of rotation of the lattice relatice to the xz plane (0-90)

% GAS7b 2D sheet dimentions
% Dist1stN=11;
% Dist2ndN=5;
% DAngle=40;
% data: Output coodinates

%% evaluate by sheet dist
pnPoints = height(data10);
boxSimu=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];


parfor ii=1:NoTrial
    if ii==1
        AvAll(:,:,ii)=(NNDfromRandomDimerBlinkNthPrintRiply3(FileName2, pnPoints,boxSimu,edges, MaxNeighbor,ObservationRate, Dist1stN,Dist2ndN,DAngle,PopulationVector, Accu, 1))';
    else
        AvAll(:,:,ii)=(NNDfromRandomDimerBlinkNthPrintRiply3(FileName2, pnPoints,boxSimu,edges, MaxNeighbor,ObservationRate, Dist1stN,Dist2ndN,DAngle,PopulationVector, Accu, 0))';       
    end
    
end

Av=mean(AvAll,3);
UpSD=Av+std(AvAll,0,3);
LwSD=Av-std(AvAll,0,3);
Av=Av';
UpSD=UpSD';
LwSD=LwSD';
end