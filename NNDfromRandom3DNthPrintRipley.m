function h3dnormV=NNDfromRandom3DNthPrintRipley(FileName2,pnPoints,boxSimu,edges, MaxNeighbor, plotflag)

simNoise = generateNoise3DtablePrint(FileName2, pnPoints,boxSimu, plotflag);
XYdataSim=[simNoise.x, simNoise.y];
ZdataSim=simNoise.z;
LdSim = height(simNoise);

% %% GPU NND
% XYdataSimGPU=gpuArray(XYdataSim);
% [IDXgpu, Dgpu]=knnsearch(XYdataSimGPU, XYdataSimGPU,'K', MaxNeighbor);
% IDX=gather(IDXgpu);
% D=gather(Dgpu);
% simNoise.IDX=IDX;
% simNoise.D=D;

[IDX, D]=knnsearch(XYdataSim, XYdataSim,'K', MaxNeighbor);
simNoise.IDX=IDX;
simNoise.D=D;
%% calculate Z distance and 3D distance

ZdifSim=zeros(pnPoints,MaxNeighbor);
D3dSim=zeros(pnPoints,MaxNeighbor);
for i=1:pnPoints 
    for iN=2:MaxNeighbor
             ZdifSim(i,iN)=ZdataSim(IDX(i,1))-ZdataSim(IDX(i,iN));
             D3dSim(i,iN)=(D(i,iN)^2+ZdifSim(i,iN)^2)^0.5;
    end
    D3dSim(i,:)=sort(D3dSim(i,:));
end
simNoise.Zdif=D3dSim;
simNoise.D3d=D3dSim;

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


end
