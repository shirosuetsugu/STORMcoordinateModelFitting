function MolSimulatedSelected=generateModelRandomDimerBlinkPrint4(FileName2,pnPoints, box, ObservationRate, Dist1stN, Dist2ndN, DAngle, PopulationVec, Accu, plotflag)
% data10: table of the original data
% pnPoints: No of data
% box=[min(data10.Xwc) max(data10.Xwc) min(data10.Ywc) max(data10.Ywc)  min(data10.Z) max(data10.Z) ];
% ObserbatioRate: % observation in the lattice (0-1)

% GAS7b 2D sheet dimentions
% Dist1stN=11;
% Dist2ndN=5;
% DAngle=40;
% data: Output coodinates

%PopulationVector ex PopulationVector=[0.52 0.27 0.145 0.062];

if nargin < 10
    plotflag=0;
end
[pathstr,name,ext] = fileparts(FileName2);

NumSimDimer=round(pnPoints/ObservationRate/2);


simNoise = generateNoise3DtablePrint(FileName2, NumSimDimer,box, 0);
XYdataSim=[simNoise.x, simNoise.y];
ZdataSim=simNoise.z;
LdSim = height(simNoise);

% generate dimer

simNoiseDimerPair=generateRandom3(NumSimDimer,Dist1stN, Accu);

simNoiseDimerPair.x=simNoiseDimerPair.x+simNoise.x;
simNoiseDimerPair.y=simNoiseDimerPair.y+simNoise.y;
simNoiseDimerPair.z=simNoiseDimerPair.z+simNoise.z;


% combine the pairs
simNoise=[simNoise; simNoiseDimerPair];

% round z because of the acurracy of the actual STORM output of z
simNoise.z=round(simNoise.z);


% random sort
DimerNum=height(simNoise);
simNoise.sort=randperm(DimerNum)';
simNoise=sortrows(simNoise,'sort');

% add the blinking probability
simNoise.prpbability=rand([DimerNum 1]);
rowsRangeP=simNoise.prpbability<=ObservationRate;
MolSimulatedSelected=simNoise(rowsRangeP,:);

% assign the  blinking numbers 6% 4 times (6.2	8.3	12.5	25)
DimerNum2=height(MolSimulatedSelected);
MolSimulatedSelected.prpbability=rand([DimerNum2 1]);
rowsRange4blinks=MolSimulatedSelected.prpbability>0 & MolSimulatedSelected.prpbability<=PopulationVec(4);
rowsRange3blinks=MolSimulatedSelected.prpbability>PopulationVec(4) & MolSimulatedSelected.prpbability<=PopulationVec(3);
rowsRange2blinks=MolSimulatedSelected.prpbability>PopulationVec(3) & MolSimulatedSelected.prpbability<=PopulationVec(2);
rowsRange1blink=MolSimulatedSelected.prpbability>PopulationVec(2) & MolSimulatedSelected.prpbability<=PopulationVec(1);
MolSimulatedSelected4blink=MolSimulatedSelected(rowsRange4blinks,:);
MolSimulatedSelected3blink=MolSimulatedSelected(rowsRange3blinks,:);
MolSimulatedSelected2blink=MolSimulatedSelected(rowsRange2blinks,:);
MolSimulatedSelected1blink=MolSimulatedSelected(rowsRange1blink,:);

% double the coordinates of 2 blinks by random 3 nm sphere
Temp=MolSimulatedSelected1blink{:,1:3};
MolSimulatedSelected1blink1=table;
MolSimulatedSelected1blink1.x=Temp(:,1);
MolSimulatedSelected1blink1.y=Temp(:,2);
MolSimulatedSelected1blink1.z=Temp(:,3);

DimerNum3=height(MolSimulatedSelected2blink);
xyz=generateRandom(DimerNum3,3);
Temp=[];Temp1=[];Temp2=[];
Temp=MolSimulatedSelected2blink{:,1:3};
centerM=(xyz{:,:})/2;
Temp1=Temp-centerM;
Temp2=Temp-centerM+xyz{:,:};
Temp3=[Temp1;Temp2];
MolSimulatedSelected2blink2=table;
MolSimulatedSelected2blink2.x=Temp3(:,1);
MolSimulatedSelected2blink2.y=Temp3(:,2);
MolSimulatedSelected2blink2.z=Temp3(:,3);

DimerNum4=height(MolSimulatedSelected3blink);
Temp=[];Temp1=[];Temp2=[];Temp3=[];
xyz=generateRandom(DimerNum4,3);
xyz2=generateRandom(DimerNum4,3);
centerM=(xyz{:,:}+xyz2{:,:})/3;
Temp=MolSimulatedSelected3blink{:,1:3};
Temp1=Temp-centerM;
Temp2=Temp-centerM+xyz{:,:};
Temp3=Temp-centerM+xyz2{:,:};
Temp4=[Temp1;Temp2;Temp3];
MolSimulatedSelected3blink3=table;
MolSimulatedSelected3blink3.x=Temp4(:,1);
MolSimulatedSelected3blink3.y=Temp4(:,2);
MolSimulatedSelected3blink3.z=Temp4(:,3);

DimerNum5=height(MolSimulatedSelected4blink);
Temp=[];Temp1=[];Temp2=[];Temp3=[];Temp4=[];
xyz=generateRandom(DimerNum5,3);
xyz2=generateRandom(DimerNum5,3);
xyz3=generateRandom(DimerNum5,3);
centerM=(xyz{:,:}+xyz2{:,:}+xyz3{:,:})/4;
Temp=MolSimulatedSelected4blink{:,1:3};
Temp1=Temp-centerM;
Temp2=Temp-centerM+xyz{:,:};
Temp3=Temp-centerM+xyz2{:,:};
Temp4=Temp-centerM+xyz3{:,:};
Temp5=[Temp1;Temp2;Temp3;Temp4];
MolSimulatedSelected4blink4=table;
MolSimulatedSelected4blink4.x=Temp5(:,1);
MolSimulatedSelected4blink4.y=Temp5(:,2);
MolSimulatedSelected4blink4.z=Temp5(:,3);

MolSimulatedSelected22=[MolSimulatedSelected4blink4;MolSimulatedSelected3blink3;MolSimulatedSelected2blink2;MolSimulatedSelected1blink1];
DimerNum6=height(MolSimulatedSelected22);
MolSimulatedSelected22.sort=randperm(DimerNum6)';
MolSimulatedSelected22=sortrows(MolSimulatedSelected22,'sort');

% add the observation probability
MolSimulatedSelected2=MolSimulatedSelected22;



% plot the results
if plotflag==1
figure
subplot(2,2,1);
line([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)]);
hold on
%scatter(simNoise.x,simNoise.y, '.');
scatter(MolSimulatedSelected2.x,MolSimulatedSelected2.y, '.');
axis equal

subplot(2,2,2);
line([box(1) box(1) box(2) box(2) box(1)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
%scatter(simNoise.x,simNoise.z, '.');
scatter(MolSimulatedSelected2.x,MolSimulatedSelected2.z, '.');
axis equal

subplot(2,2,3);
line([box(3) box(3) box(4) box(4) box(3)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
%scatter(simNoise.y,simNoise.z, '.');
scatter(MolSimulatedSelected2.y,MolSimulatedSelected2.z, '.');
axis equal

subplot(2,2,4);
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(5) box(5) box(5) box(5) box(5)], 'Color', 'b');
hold on
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(6) box(6) box(6) box(6) box(6)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(3) box(3) box(3) box(3)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(4) box(4) box(4) box(4) box(4)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
%scatter3(simNoise.x,simNoise.y,simNoise.z, '.');
scatter3(MolSimulatedSelected2.x,MolSimulatedSelected2.y,MolSimulatedSelected2.z, '.');
axis equal

FigName=[pwd filesep name '-3D-RadomDimer-NthNeighbor-Simulated-'  num2str(pnPoints) '-signals-'  num2str(Dist1stN) '-' ...
    num2str(Dist2ndN) '-' num2str( DAngle) 'lat-' num2str(ObservationRate*10) 'obs' '.pdf'];
print('-bestfit', FigName,'-dpdf','-r0');

end
end
