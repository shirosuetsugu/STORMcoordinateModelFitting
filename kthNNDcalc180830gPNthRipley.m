function [h3dnormV, h3dBinCounts]=kthNNDcalc180830gPNthRipley(FileName2, MaxNeighbor,edges, plotflag)
% Import the data from VISP
% Color by 3 more signals or 2 signals Cluster and fill
% 1:x, 2:y, 3:z, 4:dx, 5:dy, 6:dz, 7:intensity, 8: frame number
% Lxy>Lxy2, Lz>Lz2
%input parameters
if nargin < 1
    FileName2='/Users/suetsugu/Desktop/180802kfunction/mEOS4GAS7blipo-5_ZStack (30 files Z stack)4alldriftCorrTS-region4-z34320180730T192340regeionSectionRegionX20180804T154823regeionSection.txt';
end


delimiterIn = '\t';
headerlinesIn = 1;

data10=readtable(FileName2,'ReadVariableNames',true,'ReadRowNames',false ,'Delimiter','\t');

% %% GPU NND
% XYdata=[data10.Xwc, data10.Ywc];
% XYdataGPU=gpuArray(XYdata);
% Zdata=data10.Zc;
% Ld = height(data10);
% [IDXgpu, Dgpu]=knnsearch(XYdataGPU, XYdataGPU,'K', MaxNeighbor);
% IDX=gather(IDXgpu);
% D=gather(Dgpu);
% data10.IDX=IDX;
% data10.D=D;

%% CPU NND
XYdata=[data10.Xwc, data10.Ywc];
Zdata=data10.Z;
Ld = height(data10);
[IDX, D]=knnsearch(XYdata, XYdata,'K', MaxNeighbor);
data10.IDX=IDX;
data10.D=D;
%% calculate Z distance and 3D distance

Zdif=zeros(Ld,MaxNeighbor);
D3d=zeros(Ld,MaxNeighbor);
for i=1:Ld 
    for iN=2:MaxNeighbor
             Zdif(i,iN)=Zdata(IDX(i,1))-Zdata(IDX(i,iN));
             D3d(i,iN)=(D(i,iN)^2+Zdif(i,iN)^2)^0.5;
    end
    D3d(i,:)=sort(D3d(i,:));
end
data10.Zdif=Zdif;
data10.D3d=D3d;

[pathstr,name,ext] = fileparts(FileName2);
save_file_name=[pathstr filesep name '_3Ddistance.txt'];
writetable(data10,save_file_name,'WriteVariableNames',true,'WriteRowNames',false ,'Delimiter','\t');



%% volume vector
% delta=1; %bin width. Important!!
% edges = 0:delta:50;
delta=edges(2)-edges(1);
Vol1=4/3*edges.^3;
Vol0=4/3*(edges-delta).^3;
edgeL=length(edges);
Voledges=Vol1-Vol0;
Voledges=Voledges(1,1:edgeL-1);
%% histgram of 3D distance

DistColumn3D = reshape(D3d(:,2:MaxNeighbor),[],1);
h3d1=histogram(DistColumn3D,edges, 'Visible', 'off');
h3dnormV=h3d1.BinCounts./Voledges;
h3dBinCounts=h3d1.BinCounts;

if plotflag==1
    box=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];
    
    figure
    subplot(2,2,1);
    line([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)]);
    hold on
    scatter(data10.Xwc,data10.Ywc, '.');
    axis equal
    
    subplot(2,2,2);
    line([box(1) box(1) box(2) box(2) box(1)], [box(5) box(6) box(6) box(5) box(5)]);
    hold on
    scatter(data10.Xwc,data10.Z, '.');
    axis equal
    
    subplot(2,2,3);
    line([box(3) box(3) box(4) box(4) box(3)], [box(5) box(6) box(6) box(5) box(5)]);
    hold on
    scatter(data10.Ywc,data10.Z, '.');
    axis equal
    
    subplot(2,2,4);
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(5) box(5) box(5) box(5) box(5)], 'Color', 'b');
    hold on
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(6) box(6) box(6) box(6) box(6)], 'Color', 'b');
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(3) box(3) box(3) box(3)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
    plot3([box(1) box(1) box(2) box(2) box(1)], [box(4) box(4) box(4) box(4) box(4)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
    scatter3(data10.Xwc,data10.Ywc,data10.Z, '.');
    axis equal
    
    FigName=[pwd filesep name '-'  num2str(Ld) '-signals' '.pdf'];
    print('-bestfit', FigName,'-dpdf','-r0');
    
end


end