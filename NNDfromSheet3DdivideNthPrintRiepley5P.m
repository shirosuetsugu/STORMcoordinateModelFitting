function h3dnormV=NNDfromSheet3DdivideNthPrintRiepley5P(FileName2, ...
    pnPoints,box,edges, MaxNeighbor,ObservationRate, Dist1stN,Dist2ndN,DAngle, tiltAngleX, tiltAngleZ, divNum,PopulationVector, Accu,PercentSheet, plotflag)

simSheet1=generateModelLatticeParallelDividePrint6(FileName2, round(pnPoints*PercentSheet), box, ObservationRate, Dist1stN,Dist2ndN,DAngle, tiltAngleX, tiltAngleZ,divNum,PopulationVector,Accu,  plotflag);
simSheet2=generateModelRandomDimerBlinkPrint4(FileName2, round(pnPoints*(1-PercentSheet)), box, ObservationRate, Dist1stN,Dist2ndN,DAngle,PopulationVector,Accu, plotflag);
simSheet=table;
simSheet.x=[simSheet1.x;simSheet2.x];
simSheet.y=[simSheet1.y;simSheet2.y];
simSheet.z=[simSheet1.z;simSheet2.z];

XYdataSim=[simSheet.x, simSheet.y];
ZdataSim=simSheet.z;
LdSim = height(simSheet);

% %% GPU NND
% XYdataSimGPU=gpuArray(XYdataSim);
% [IDXgpu, Dgpu]=knnsearch(XYdataSimGPU, XYdataSimGPU,'K', MaxNeighbor);
% IDX=gather(IDXgpu);
% D=gather(Dgpu);
% simNoise.IDX=IDX;
% simNoise.D=D;

[IDX, D]=knnsearch(XYdataSim, XYdataSim,'K', MaxNeighbor);
simSheet.IDX=IDX;
simSheet.D=D;
%% calculate Z distance and 3D distance

ZdifSim=zeros(LdSim,MaxNeighbor);
D3dSim=zeros(LdSim,MaxNeighbor);
for i=1:LdSim 
    for iN=2:MaxNeighbor
             ZdifSim(i,iN)=ZdataSim(IDX(i,1))-ZdataSim(IDX(i,iN));
             D3dSim(i,iN)=(D(i,iN)^2+ZdifSim(i,iN)^2)^0.5;
    end
    D3dSim(i,:)=sort(D3dSim(i,:));
end
simSheet.Zdif=D3dSim;
simSheet.D3d=D3dSim;

%% volume vector
delta=edges(2)-edges(1);
Vol1=4/3*edges.^3;
Vol0=4/3*(edges-delta).^3;
edgeL=length(edges);
Voledges=Vol1-Vol0;
Voledges=Voledges(1,1:edgeL-1);
%% histgram of 3D distance
% plot histgrams

DistColumn3DSim = reshape(D3dSim(:,2:MaxNeighbor),[],1);
h3d1=histogram(DistColumn3DSim,edges, 'Visible', 'off');
h3dnormV=h3d1.BinCounts./Voledges;
h3dBinCounts=h3d1.BinCounts;

if plotflag==1
    figure
    subplot(2,2,1);
    line([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)]);
    hold on
    scatter(simSheet.x,simSheet.y, '.');
    axis equal
    
    subplot(2,2,2);
    line([box(1) box(1) box(2) box(2) box(1)], [box(5) box(6) box(6) box(5) box(5)]);
    hold on
    scatter(simSheet.x,simSheet.z, '.');
    axis equal
    
    subplot(2,2,3);
    line([box(3) box(3) box(4) box(4) box(3)], [box(5) box(6) box(6) box(5) box(5)]);
    hold on
    scatter(simSheet.y,simSheet.z, '.');
    axis equal
    
    subplot(2,2,4);
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(5) box(5) box(5) box(5) box(5)], 'Color', 'b');
    hold on
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(6) box(6) box(6) box(6) box(6)], 'Color', 'b');
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(3) box(3) box(3) box(3)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(4) box(4) box(4) box(4) box(4)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
    scatter3(simSheet.x,simSheet.y,simSheet.z, '.');
    axis equal
    [pathstr,name,ext] = fileparts(FileName2);
    FigName=[pwd filesep name '-3D-NthNeighbor-Simulated-MixPSheet-'  num2str(pnPoints) '-signals-'  num2str(Dist1stN) '-' ...
        num2str(Dist2ndN) '-' num2str( DAngle) 'lat-' num2str( tiltAngleX) '-' num2str( tiltAngleZ) 'tilt-' num2str(ObservationRate*10) 'obs' '.pdf'];
    print('-bestfit', FigName,'-dpdf','-r0');
    
end

end
