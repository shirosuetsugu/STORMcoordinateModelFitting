function EvaluateDimerDist180923RipleyLParallelModel(FileName2, MaxNeighbor, Dist1stN,Dist2ndN,DAngle, NoTrial, PopulationVec, Para1, edges, Accu, PercentSheet)

%GAS7 lattice
% Dist1stN=11;
% Dist2ndN=5;
% DAngle=40;
% 
% MaxNeighbor=5; %+1
% 
% NoTrial=5; %for fig 20
% 
% %%tiltAngleX, Z, Observation Rate->Para1
% ObservationRate=0.5;
% tiltAngleX=40;
% tiltAngleZ=40;

%% Parameters
% Para1=table;
% ParaNum=3;
% Para1.tiltAngleX=zeros(ParaNum,1);
% Para1.tiltAngleZ=zeros(ParaNum,1);
% Para1.ObservationRate=[0.075 0.1 0.2]';
% Para1.DivNum=3;
ParaNum=height(Para1);

%%
ParaSize=height(Para1);

[pathstr,name,ext] = fileparts(FileName2);


NewDirName=[pathstr filesep name 'RandomDimerSimulatedBlinkRipley' num2str(MaxNeighbor) ];
status = mkdir(NewDirName);
cd(NewDirName);

%% data
% delimiterIn = '\t';
% headerlinesIn = 1;

data10=readtable(FileName2,'ReadVariableNames',true,'ReadRowNames',false ,'Delimiter','\t');
pnPoints = height(data10);
box=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];
VolBox=(box(2)-box(1))*(box(4)-box(3))*(box(6)-box(5));
MolConcSignal=pnPoints/6.02e23/VolBox/1e-24;

% data values
[h3dnormV, h3dBinCounts]=kthNNDcalc180830gPNthRipley(FileName2, MaxNeighbor,edges,1);

% simulated random dimers
for kk=1:ParaNum
    fprintf('%s\r', num2str(kk));
    [SimAv, SimUpSD, SimLwSD, Av, UpSD, LwSD]=kthNNDcalcSimulated180918RandomDimerBlinkNthRipley(FileName2, ...
        MaxNeighbor,Para1.ObservationRate(kk), Dist1stN,Dist2ndN,DAngle, NoTrial, PopulationVec, edges, Accu);
    SimAvAll(:,:,kk)=SimAv;
    SimUpSDAll(:,:,kk)=SimUpSD;
    SimLwSDAll(:,:,kk)=SimLwSD;
    AvAll(:,:,kk)=Av;
    UpSDAll(:,:,kk)=UpSD;
    LwSDAll(:,:,kk)=LwSD;
    edgesAll(:,:,kk)=edges;
    close all
end

% simulated sheets 
for kk=1:ParaNum
    fprintf('%s\r', num2str(kk));
    [SimAvSheet, SimUpSDSheet, SimLwSDSheet, AvSheet, UpSDSheet, LwSDSheet]=kthNNDcalcSimulated180922gPNthRiepleyP(FileName2, ...
        MaxNeighbor, Para1.ObservationRate(kk), Dist1stN,Dist2ndN,DAngle, Para1.tiltAngleX(kk), Para1.tiltAngleZ(kk), NoTrial, Para1.DivNum(kk),  PopulationVec, edges, Accu,PercentSheet);
    SimAvSheetAll(:,:,kk)=SimAvSheet;
    SimUpSDSheetAll(:,:,kk)=SimUpSDSheet;
    SimLwSDSheetAll(:,:,kk)=SimLwSDSheet;
    AvSheetAll(:,:,kk)=AvSheet;
    UpSDSheetAll(:,:,kk)=UpSDSheet;
    LwSDSheetAll(:,:,kk)=LwSDSheet;
    edgesAll(:,:,kk)=edges;
    close all
end


FileNameSave=[pwd filesep name  'Simulated.mat'];
save(FileNameSave, 'SimAvAll','AvAll','-v7.3');

%% average by grouping of Observation rate (given by Para)
edges=edgesAll(:,:,1);
for jj=1:1 %this value depends on Para file
    
    edgesL=length(edges);
    edgesDL=length(edges);
%% Figure 1
    figure
    hold on
    for kk=1:ParaNum
        c = colorcode(kk+1);
      %  scatter(edges(1,4:edgeL-1)',SimAvAll(:,2:end,kk)./AvAll(:,2:end,kk),'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
        plot(edges(1,2:edgesL-1)',SimAvAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-', 'Color' ,c ,'LineWidth',2);
        plot(edges(1,2:edgesL-1)',SimUpSDAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',1);
        plot(edges(1,2:edgesL-1)',SimLwSDAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',1);
        plot(edges(1,2:edgesL-1)',AvAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle',':', 'Color' ,c ,'LineWidth',0.5);
        plot(edges(1,2:edgesL-1)',UpSDAll(:,2:end,kk)./AvAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.5);
        plot(edges(1,2:edgesL-1)',LwSDAll(:,2:end,kk)./AvAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.5);
    end
    
    for kk=1:ParaNum
        c = colorcode(kk+1);
        %  scatter(edges(1,4:edgeL-1)',SimAvSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk),'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
        plot(edges(1,2:edgesL-1)',SimAvSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk), 'LineStyle', '-', 'Color' ,c ,'LineWidth',1.5);
        plot(edges(1,2:edgesL-1)',SimUpSDSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',0.75);
        plot(edges(1,2:edgesL-1)',SimLwSDSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',0.75);
        plot(edges(1,2:edgesL-1)',AvSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk), 'LineStyle',':', 'Color' ,c ,'LineWidth',0.25);
        plot(edges(1,2:edgesL-1)',UpSDSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.25);
        plot(edges(1,2:edgesL-1)',LwSDSheetAll(:,2:end,kk)./AvSheetAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.25);
    end
    
    %  plot data ratio
    for kk=1:ParaNum
        h3dnormVRatioAll(:,:,kk)=h3dnormV(:,2:end)./AvAll(:,2:end,kk);
                
    end
    h3dnormVplotAll=mean(h3dnormVRatioAll,3);
    h3dnormVplotAllUpSD=h3dnormVplotAll+std(h3dnormVRatioAll,0,3);
    h3dnormVplotAllLwSD=h3dnormVplotAll-std(h3dnormVRatioAll,0,3);
    c = colorcode(2);
    hold on
    %         scatter(edges(1,2:edgesL-1),h3dnormV(:,2:end)./AvAll(:,2:end,kk),'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
    plot(edges(1,2:edgesL-1),h3dnormVplotAll(:,:), 'LineStyle', '-', 'Color' ,c ,'LineWidth',3 );
    plot(edges(1,2:edgesL-1),h3dnormVplotAllUpSD(:,:), 'LineStyle', '-.', 'Color' ,c ,'LineWidth',3 );
    plot(edges(1,2:edgesL-1),h3dnormVplotAllLwSD(:,:), 'LineStyle', '-.', 'Color' ,c ,'LineWidth',3 );
    
    Title1=[name ];
    Title11=' 3D-NthNeighbor-Simulated-ratio-ParallelSheet- ';
    Title2=[  num2str(MaxNeighbor) 'MaxNeighbors-'  num2str(pnPoints) 'signals-'  num2str(Dist1stN) 'nm-' ...
        num2str(Dist2ndN) 'nm-' num2str( DAngle) 'Angle-TiltX' num2str( Para1.tiltAngleX(jj)) 'TiltZ-'  num2str(Para1.tiltAngleZ(jj))];
    Title3=['Accuracy ' num2str(Accu) '-' num2str(PercentSheet*100) '%sheet, Signalconc ' num2str(MolConcSignal) ' mol/l' ];
    Title4=['PopulationVector ' num2str(PopulationVec)];
    Title4=Title4(~ isspace (Title4));
    Title5=num2str(Para1.ObservationRate');
    Title5=Title5(~ isspace (Title5));
    title({Title1, Title11, Title2, Title3, Title4 , Title5})
    
      
    %     ylim([-10 300]);
    FigName=[pwd filesep name '-3D-data-RadomDimer-SimulatedParallelSheet-' num2str(pnPoints) 'signals' num2str(Dist1stN) '-' num2str(Dist2ndN) '-'  num2str( DAngle) '-'  num2str( Para1.tiltAngleX(jj)) '-' num2str(Para1.tiltAngleZ(jj)) '-' num2str(jj) '-' Title4 Title5 '.pdf'];
    print('-bestfit', FigName,'-dpdf','-r0');
    
    %% Figure 2
    figure
    hold on
    for kk=1:ParaNum
        c = colorcode(kk+1);
      %  scatter(edges(1,4:edgeL-1)',SimAvAll(:,2:end,kk)./AvAll(:,2:end,kk),'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
        plot(edges(1,2:edgesL-1)',SimAvAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-', 'Color' ,c ,'LineWidth',2);
        plot(edges(1,2:edgesL-1)',SimUpSDAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',1);
        plot(edges(1,2:edgesL-1)',SimLwSDAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle', '-.', 'Color' ,c  ,'LineWidth',1);
        plot(edges(1,2:edgesL-1)',AvAll(:,2:end,kk)./AvAll(:,2:end,kk), 'LineStyle',':', 'Color' ,c ,'LineWidth',0.5);
        plot(edges(1,2:edgesL-1)',UpSDAll(:,2:end,kk)./AvAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.5);
        plot(edges(1,2:edgesL-1)',LwSDAll(:,2:end,kk)./AvAll(:,2:end,kk),'LineStyle','-.', 'Color' ,c,'LineWidth',0.5);
    end
    
   
    %  plot data ratio
    for kk=1:ParaNum
        h3dnormVRatioAll(:,:,kk)=h3dnormV(:,2:end)./AvAll(:,2:end,kk);
                
    end
    h3dnormVplotAll=mean(h3dnormVRatioAll,3);
    h3dnormVplotAllUpSD=h3dnormVplotAll+std(h3dnormVRatioAll,0,3);
    h3dnormVplotAllLwSD=h3dnormVplotAll-std(h3dnormVRatioAll,0,3);
    c = colorcode(2);
    hold on
    %         scatter(edges(1,2:edgesL-1),h3dnormV(:,2:end)./AvAll(:,2:end,kk),'MarkerEdgeColor', c , 'MarkerFaceColor', 'none' );
    plot(edges(1,2:edgesL-1),h3dnormVplotAll(:,:), 'LineStyle', '-', 'Color' ,c ,'LineWidth',3 );
    plot(edges(1,2:edgesL-1),h3dnormVplotAllUpSD(:,:), 'LineStyle', '-.', 'Color' ,c ,'LineWidth',3 );
    plot(edges(1,2:edgesL-1),h3dnormVplotAllLwSD(:,:), 'LineStyle', '-.', 'Color' ,c ,'LineWidth',3 );

    Title1=[name ];
    Title11=' 3D-NthNeighbor-Simulated-ratio-randomdimer- ';
    Title2=[  num2str(MaxNeighbor) 'MaxNeighbors-'  num2str(pnPoints) 'signals-'  num2str(Dist1stN) 'nm-' ...
        num2str(Dist2ndN) 'nm-' num2str( DAngle) 'Angle-TiltX' num2str( Para1.tiltAngleX(jj)) 'TiltZ-'  num2str(Para1.tiltAngleZ(jj))];
    Title3=[ 'Signalconc ' num2str(MolConcSignal) ' mol/l' ];
    Title4=['PopulationVector ' num2str(PopulationVec)];
    Title4=Title4(~ isspace (Title4));
    Title5=num2str(Para1.ObservationRate');
    Title5=Title5(~ isspace (Title5));
    title({Title1, Title11, Title2, Title3, Title4 , Title5})
    
      
    %     ylim([-10 300]);
    FigName=[pwd filesep name '-3D-data-RadomDimer-' num2str(pnPoints) 'signals' num2str(Dist1stN) '-' num2str(Dist2ndN) '-'  num2str( DAngle) '-' num2str( Para1.tiltAngleX(jj)) '-' num2str(Para1.tiltAngleZ(jj)) '-' num2str(jj) '-' Title4 Title5 '.pdf'];
    print('-bestfit', FigName,'-dpdf','-r0');
    
    %% plot histgrams
    formatOut= 30;
    c=datestr(now, formatOut);
    
    figure
    hold on
    c = colorcode(2);
    plot(edges(1,1:edgesDL-1),h3dBinCounts,'Color' ,c);
    
    Title1=[name ' 3DHistDensity '  num2str(MaxNeighbor) ];
    title(Title1)
    FigName=[pwd filesep name '-3DHist-total'  num2str(pnPoints) 'signals'  '.pdf'];
    print('-bestfit', FigName,'-dpdf','-r0');
    
    
    %% plot with random probabilities
    %
    % figure
    % hold on
    %
    %     c = colorcode(2);
    %     scatter(edges(1,2:edgesL-1),h3dnormV,'MarkerEdgeColor' ,c, 'MarkerFaceColor', 'none' );
    %     plot(edges(1,2:edgesL-1),h3dnormV, 'LineStyle', '-', 'Color' ,c );
    %     plot(edges(1,2:edgesL-1),Av, 'LineStyle',':', 'Color' ,c );
    %     plot(edges(1,2:edgesL-1),UpSD,'LineStyle','-.', 'Color' ,c );
    %     plot(edges(1,2:edgesL-1),LwSD,'LineStyle','-.', 'Color' ,c );
    %
    % ax = gca;
    % ax.YScale = 'log';
    %
    % Title1=[name ' 3D-NthNeighbor '  num2str(MaxNeighbor) ];
    % title(Title1)
    % FigName=[pwd filesep name '-3D-NthNeighbor-'  num2str(pnPoints) 'signals'  '.pdf'];
    % print('-bestfit', FigName,'-dpdf','-r0');
    
end

%close all

