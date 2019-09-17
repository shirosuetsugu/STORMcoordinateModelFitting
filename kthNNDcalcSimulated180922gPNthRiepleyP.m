function [SimAv, SimUpSD, SimLwSD, Av, UpSD, LwSD]=kthNNDcalcSimulated180922gPNthRiepleyP(FileName2, ...
    MaxNeighbor, ObservationRate, Dist1stN,Dist2ndN,DAngle, tiltAngleX, tiltAngleZ, NoTrial, divNum,  PopulationVec, edges,Accu, PercentSheet)
% ObservationRate=0.5;
% tiltAngleX=40;
% tiltAngleZ=40;
% 
% Dist1stN=11;
% Dist2ndN=5;
% DAngle=40;

% MaxNeighbor=5; %+1



% FileName2='/Users/suetsugu/Desktop/180802kfunction/mEOS4GAS7blipo-5_ZStack (30 files Z stack)4alldriftCorrTS-region4-z34320180730T192340regeionSectionRegionX20180804T154823regeionSection.txt';
delimiterIn = '\t';
headerlinesIn = 1;
[pathstr,name,ext] = fileparts(FileName2);

data10=readtable(FileName2,'ReadVariableNames',true,'ReadRowNames',false ,'Delimiter','\t');
Ld = height(data10);

%% lattice generation
 % No of simulations 

pnPoints = height(data10);
box=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];
% delta=1; %bin width. Important!!
% edges = 0:delta:50;]
delta=edges(2)-edges(1);
Vol1=4/3*edges.^3;
Vol0=4/3*(edges-delta).^3;
edgeL=length(edges);

[SimAv, SimUpSD, SimLwSD]=kthNNDcalcSheetNth180922RieplyP(FileName2, data10, MaxNeighbor, edges, NoTrial,ObservationRate, Dist1stN,Dist2ndN,DAngle, tiltAngleX, tiltAngleZ, divNum,PopulationVec, Accu, PercentSheet);

%% random probabilities
[Av, UpSD, LwSD]=kthNNDcalcRandomNth180830Ripley(FileName2, data10, MaxNeighbor, edges, NoTrial);

%%




figure
hold on
c = colorcode(2);
scatter(edges(1,1:edgeL-1),SimAv./Av,'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
%     TempPlot1(~isfinite(h3dnormV./Av))=0;
plot(edges(1,1:edgeL-1),SimAv./Av, 'LineStyle', '-', 'Color' ,c  );
plot(edges(1,1:edgeL-1),SimUpSD./Av, 'LineStyle', '-.', 'Color' ,c  );
plot(edges(1,1:edgeL-1),SimLwSD./Av, 'LineStyle', '-.', 'Color' ,c  );
plot(edges(1,1:edgeL-1),Av./Av, 'LineStyle',':', 'Color' ,c );
plot(edges(1,1:edgeL-1),UpSD./Av,'LineStyle','-.', 'Color' ,c);
plot(edges(1,1:edgeL-1),LwSD./Av,'LineStyle','-.', 'Color' ,c);



Title1=[name ' 3D-NthNeighbor-Simulated-Sheet- '];
Title2=[  num2str(MaxNeighbor) 'MN-'  num2str(Ld) 'signals-' num2str(ObservationRate) 'ObsR-' num2str(Dist1stN) 'nm-' ...
    num2str(Dist2ndN) 'nm-' num2str( DAngle) 'degTilt-' num2str( tiltAngleX) 'degXrot-' num2str( tiltAngleZ) 'degZrot'];
title({Title1, Title2})

FigName=[pwd filesep name '-3D-NthNeighbor-Simulated-Sheet-'  num2str(Ld) '-signals-' num2str(ObservationRate) '-' num2str(Dist1stN) '-' ...
    num2str(Dist2ndN) '-' num2str( DAngle) '-' num2str( tiltAngleX) '-' num2str( tiltAngleZ) '.pdf'];
print('-bestfit', FigName,'-dpdf','-r0');
