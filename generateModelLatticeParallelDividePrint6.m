function MolSimulatedSelected=generateModelLatticeParallelDividePrint6(FileName2,pnPoints, box, ObservationRate, Dist1stN, Dist2ndN, DAngle, tiltAngleX, tiltAngleZ,divNum,PopulationVec, Accu, plotflag)
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

if nargin < 11
    plotflag=0;
end
[pathstr,name,ext] = fileparts(FileName2);

NumSimData=round(pnPoints/ObservationRate);
centerMass=[mean(box(1:2)),mean(box(3:4)),mean(box(5:6))];


MolArea=Dist1stN*Dist2ndN*sin(DAngle/180*pi());
SheetArea=MolArea*NumSimData;
SideLength=sqrt(SheetArea);

i1stMax=round(SideLength/Dist1stN);
i2ndMax=round(SideLength/(Dist2ndN*sin(DAngle/180*pi())));
TempNumData=round(max([box(2)-box(1),box(4)-box(3)])^2*2/MolArea);
TempSideLength=sqrt(MolArea*TempNumData);

MolSimulated=table;
MolSimulated.x=zeros(TempNumData,1);
MolSimulated.y=zeros(TempNumData,1);
MolSimulated.z=zeros(TempNumData,1);
BackFactor=round(Dist1stN/(Dist2ndN*sin(DAngle/180*pi())));

for ii=1:TempNumData
    if ii==1
        MolSimulated.x(ii)=min(box(1:2));
        MolSimulated.y(ii)=min(box(3:4));
        MolSimulated.z(ii)=mean(box(5:6));
        jj=ii ;
    elseif ii>1 && MolSimulated.x(ii-1)<box(1)+(box(2)-box(1))/2+SideLength/2
        MolSimulated.x(ii)=MolSimulated.x(ii-1)+Dist1stN;
        MolSimulated.y(ii)=MolSimulated.y(ii-1);
        MolSimulated.z(ii)=MolSimulated.z(ii-1);
    elseif ii>1 && MolSimulated.x(ii-1)>=box(1)+(box(2)-box(1))/2+SideLength/2 && ~mod(jj,BackFactor)==0
        MolSimulated.x(ii)=MolSimulated.x(jj)+Dist2ndN*cos(DAngle/180*pi());
        MolSimulated.y(ii)=MolSimulated.y(jj)+Dist2ndN*sin(DAngle/180*pi());
        MolSimulated.z(ii)=MolSimulated.z(jj);
        jj=ii ;
    elseif ii>1 && MolSimulated.x(ii-1)>=box(1)+(box(2)-box(1))/2+SideLength/2 && mod(jj,BackFactor)==0
        MolSimulated.x(ii)=MolSimulated.x(jj)+Dist2ndN*cos(DAngle/180*pi())-Dist1stN;
        MolSimulated.y(ii)=MolSimulated.y(jj)+Dist2ndN*sin(DAngle/180*pi());
        MolSimulated.z(ii)=MolSimulated.z(jj);
        jj=ii ;

    end
end
% center of mass
TempCoord.x=MolSimulated.x-mean(MolSimulated.x);
TempCoord.y=MolSimulated.y-mean(MolSimulated.y);
TempCoord.z=MolSimulated.z-mean(MolSimulated.z);


% rotate at Z axis
tZ=tiltAngleZ/180*pi();
Rz = [cos(tZ) -sin(tZ) 0; sin(tZ) cos(tZ) 0; 0 0 1];
% rotate at X axis
tX=tiltAngleX/180*pi();
Rx = [1 0 0; 0 cos(tX) -sin(tX); 0 sin(tX) cos(tX)];
xyz=[TempCoord.x';TempCoord.y';TempCoord.z'];
xyzRx = Rz * (Rx * xyz);

MolSimulated.x=xyzRx(1,:)'+centerMass(1);
MolSimulated.y=xyzRx(2,:)'+centerMass(2);
MolSimulated.z=xyzRx(3,:)'+centerMass(3);


%move the region outside the range (only at Z direction)
BoxThickness=box(6)-box(5);
rowsRangeTemp1=[];rowsRangeTemp2=[];
while min(MolSimulated.z)<box(5)
    rowsRangeTemp1=MolSimulated.z<box(5);
    rowsRangeTemp2=MolSimulated.z>=box(5);
    Temp1=MolSimulated(rowsRangeTemp1,:);
    Temp1.z=Temp1.z+BoxThickness;
    Temp1.x=Temp1.x-2*Dist1stN; % avoid overlap after moving
    Temp1.y=Temp1.y-2*Dist2ndN;% avoid overlap after moving
    Temp2=MolSimulated(rowsRangeTemp2,:);
    MolSimulated=[Temp1;Temp2];
end
rowsRangeTemp1=[];rowsRangeTemp2=[];
while max(MolSimulated.z)>box(6)
    rowsRangeTemp1=MolSimulated.z>box(6);
    rowsRangeTemp2=MolSimulated.z<=box(6);
    Temp1=MolSimulated(rowsRangeTemp1,:);
    Temp1.z=Temp1.z-BoxThickness;
    Temp1.x=Temp1.x+2*Dist1stN; % avoid overlap after moving
    Temp1.y=Temp1.y+2*Dist2ndN; % avoid overlap after moving
    Temp2=MolSimulated(rowsRangeTemp2,:);
    MolSimulated=[Temp1;Temp2];
end

% Take required no of data
MolSimulated.number=(1:height(MolSimulated))';
rowsRange12=MolSimulated.x>=box(1) &  MolSimulated.x<=box(2) & MolSimulated.y>=box(3) & MolSimulated.y<=box(4);
MolSimulated1=MolSimulated(rowsRange12,:);
if NumSimData<height(MolSimulated)
    MolSimulated1=MolSimulated(1:NumSimData,:);
end

% move to center
MolSimulated11=table;
MolSimulated11.x=MolSimulated1.x-mean(MolSimulated1.x)+centerMass(1);
MolSimulated11.y=MolSimulated1.y-mean(MolSimulated1.y)+centerMass(2);
MolSimulated11.z=MolSimulated1.z-mean(MolSimulated1.z)+centerMass(3);

% round z because of the acurracy of the actual STORM output of z
MolSimulated11.z=round(MolSimulated1.z);





% dived by divNum
boxSimu=[min(MolSimulated11.x)-1 max(MolSimulated11.x)+1 min(MolSimulated11.y)-1 max(MolSimulated11.y)+1  min(MolSimulated11.z)-1 max(MolSimulated11.z)+1 ];
XsimDivPoint=boxSimu(1):(boxSimu(2)-boxSimu(1))/divNum:boxSimu(2);
YsimDivPoint=boxSimu(3):(boxSimu(4)-boxSimu(3))/divNum:boxSimu(4);
ZsimDivPoint=boxSimu(5):(boxSimu(6)-boxSimu(5))/1:boxSimu(6);
XdivPoint=box(1):(box(2)-box(1))/divNum:box(2);
YdivPoint=box(3):(box(4)-box(3))/divNum:box(4);
ZdivPoint=box(5):(box(6)-box(5))/1:box(6);

totalDiv=divNum^2;
XYZsimdivPoints=zeros(totalDiv,6);
XYZdivPoints=zeros(totalDiv,6);

for ll=1:totalDiv
    XYZsimdivPoints(ll,1)=XsimDivPoint(fix((ll-1)/divNum)+1);
    XYZsimdivPoints(ll,2)=XsimDivPoint(fix((ll-1)/divNum)+1+1);
    XYZsimdivPoints(ll,3)=YsimDivPoint(ll-divNum*fix((ll-1)/divNum));
    XYZsimdivPoints(ll,4)=YsimDivPoint(ll+1-divNum*fix((ll-1)/divNum));
    XYZsimdivPoints(ll,5)=ZsimDivPoint(1);
    XYZsimdivPoints(ll,6)=ZsimDivPoint(2);
    
    XYZdivPoints(ll,1)=XdivPoint(fix((ll-1)/divNum)+1);
    XYZdivPoints(ll,2)=XdivPoint(fix((ll-1)/divNum)+1+1);
    XYZdivPoints(ll,3)=YdivPoint(ll-divNum*fix((ll-1)/divNum));
    XYZdivPoints(ll,4)=YdivPoint(ll+1-divNum*fix((ll-1)/divNum));
    XYZdivPoints(ll,5)=ZdivPoint(1);
    XYZdivPoints(ll,6)=ZdivPoint(2);
end

MolSimulatedDiv=table;
for mm=1:totalDiv
    rowsRangeTemp11=MolSimulated11.x>=XYZsimdivPoints(mm,1) & MolSimulated11.x<XYZsimdivPoints(mm,2) ...
        & MolSimulated11.y>=XYZsimdivPoints(mm,3) & MolSimulated11.y<XYZsimdivPoints(mm,4) ...
        & MolSimulated11.z>=XYZsimdivPoints(mm,5) & MolSimulated11.z<XYZsimdivPoints(mm,6);
    Temp21=MolSimulated11(rowsRangeTemp11,:);
    
    Temp21.x=Temp21.x-XYZsimdivPoints(mm,1)+XYZdivPoints(mm,1);
    Temp21.y=Temp21.y-XYZsimdivPoints(mm,3)+XYZdivPoints(mm,3);
    Temp21.z=Temp21.z-XYZsimdivPoints(mm,5)+XYZdivPoints(mm,5);
    MolSimulatedDiv=[MolSimulatedDiv;Temp21];
end

MolSimulated=MolSimulatedDiv;
MolSimulated.x=MolSimulated.x-mean(MolSimulated.x)+centerMass(1);
MolSimulated.y=MolSimulated.y-mean(MolSimulated.y)+centerMass(2);
MolSimulated.z=MolSimulated.z-mean(MolSimulated.z)+centerMass(3);


% random sort
MolSimulated.sort=randperm(height(MolSimulated))';
MolSimulated=sortrows(MolSimulated,'sort');

% add the probability
MolSimulated.prpbability=rand([height(MolSimulated) 1]);
rowsRangeP=MolSimulated.prpbability<=ObservationRate;
MolSimulatedSelected=MolSimulated(rowsRangeP,:);


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

% move to the center
xMove=mean(MolSimulatedSelected22.x)-centerMass(1);
yMove=mean(MolSimulatedSelected22.y)-centerMass(2);
zMove=mean(MolSimulatedSelected22.z)-centerMass(3);

MolSimulatedSelected22.x=MolSimulatedSelected22.x-xMove;
MolSimulatedSelected22.y=MolSimulatedSelected22.y-yMove;
MolSimulatedSelected22.z=MolSimulatedSelected22.z-zMove;


% add the observation probability
MolSimulatedSelected2=MolSimulatedSelected22;


% add the spatial rondomness
simNoiseAdd=generateRandom3(height(MolSimulatedSelected2),0, Accu);

MolSimulatedSelected2.x=MolSimulatedSelected2.x+simNoiseAdd.x;
MolSimulatedSelected2.y=MolSimulatedSelected2.y+simNoiseAdd.y;
MolSimulatedSelected2.z=MolSimulatedSelected2.z+simNoiseAdd.z;

%%%%%%
% plot the results
if plotflag==1
figure
subplot(2,2,1);
line([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)]);
hold on
scatter(MolSimulated.x,MolSimulated.y, '.');
scatter(MolSimulatedSelected2.x,MolSimulatedSelected2.y, '.');
axis equal

subplot(2,2,2);
line([box(1) box(1) box(2) box(2) box(1)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
scatter(MolSimulated.x,MolSimulated.z, '.');
scatter(MolSimulatedSelected2.x,MolSimulatedSelected2.z, '.');
axis equal

subplot(2,2,3);
line([box(3) box(3) box(4) box(4) box(3)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
scatter(MolSimulated.y,MolSimulated.z, '.');
scatter(MolSimulatedSelected2.y,MolSimulatedSelected2.z, '.');
axis equal

subplot(2,2,4);
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(5) box(5) box(5) box(5) box(5)], 'Color', 'b');
hold on
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(6) box(6) box(6) box(6) box(6)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(3) box(3) box(3) box(3)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(4) box(4) box(4) box(4) box(4)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
scatter3(MolSimulated.x,MolSimulated.y,MolSimulated.z, '.');
scatter3(MolSimulatedSelected2.x,MolSimulatedSelected2.y,MolSimulatedSelected2.z, '.');
axis equal

FigName=[pwd filesep name '-3D-NthNeighbor-Simulated-Sheet-'  num2str(pnPoints) '-signals-'  num2str(Dist1stN) '-' ...
    num2str(Dist2ndN) '-' num2str( DAngle) 'lat-' num2str( tiltAngleX) '-' num2str( tiltAngleZ) 'tilt-' num2str(ObservationRate*10) 'obs' '.pdf'];
print('-bestfit', FigName,'-dpdf','-r0');

end
end

