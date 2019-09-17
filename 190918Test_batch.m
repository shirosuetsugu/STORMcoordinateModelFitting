%test for batch job on the server
DirName='/Users/suetsugu/Desktop/180924kfunctionFinal';
%DirName='/work/suetsugu/190902kfunctionFinal';
FileNameT{1}=[ DirName filesep '160914_5_2_C2_002_list-2016-12-27-09-16-07_S01corrMinPhoto100driftCorrTS-region2-z26820180730T191641regeionSectionregeionSection-28.txt' ]; % Fig 4i
FileNameT{2}=[ DirName filesep 'mEOS4GAS7blipo-5_ZStack (30 files Z stack)4alldriftCorrTS-region4-z34320180730T192340regeionSectionRegionX20180804T154823regeionSection.txt' ]; % Fig 4g

FileNameT{3}=[ DirName filesep 'mEOS4GAS7blipo-5_ZStack (30 files Z stack)4alldriftCorrTS-region-blank1regeionSectionregeionSection-11.txt' ]; % Fig 4f
FileNameT{4}=[ DirName filesep '160914_5_2_C2_002_list-2016-12-27-09-16-07_S01corrMinPhoto100driftCorrTS-region2-z4320180904T140415regeionSectionregeionSection-33.txt' ]; % Fig 4h


MaxNeighbor=10; %+1

NoTrial=20; %for fig 20

%% population of the molecules with 1, 2, 3, 4 signals
PopulationVector{1}=[0.52 0.27 0.145 0.062];
PopulationVector{2}=[0.75 0.25 0 0];
PopulationVector{3}=[0.9 0.1 0 0];
PopulationVector{4}=[0.73 0.2 0.07 0];
PopulationVector{5}=[1 0 0 0];
PopulationVector{6}=[0.50 0.27 0.145 0.1];
PopulationVector{7}=[0.3 0.3 0.2 0.1];

%% define the angle of observation and observation rate
Para1=table;
ParaNum=3;
Para1.tiltAngleX=[0 30 60 ]';
Para1.tiltAngleZ=zeros(ParaNum,1);
Para1.ObservationRate=[0.05 0.10 0.15]';
Para1.DivNum=1*ones(ParaNum,1);

delta=5; %bin width. Important!!
edges = 0:delta:50;
Accu=20; %accuracy 20 nm
PercentSheet=0.4;

%% define the GAS7 lattice size
Dist1stN=11;
Dist2ndN=5;
DAngle=40;


%% Simulate and write graph
for  ii3=1:4
    [pathstr,name,ext] = fileparts(FileNameT{ii3});
    fprintf('%s\r', num2str(ii3));
    fprintf('%s\r', name);
    
    for tt=5:5
        EvaluateDimerDist180923RipleyLParallelModel(FileNameT{ii3}, MaxNeighbor, Dist1stN,Dist2ndN,DAngle, NoTrial, PopulationVector{tt}, Para1, edges, Accu, PercentSheet);
    end
end

for  ii3=1:4
    [pathstr,name,ext] = fileparts(FileNameT{ii3});
    fprintf('%s\r', num2str(ii3));
    fprintf('%s\r', name);
    
    for tt=5:5 % FFO organization is random
        EvaluateDimerDist180923RipleyLAntiParallelModel(FileNameT{ii3}, MaxNeighbor, Dist1stN,Dist2ndN,DAngle, NoTrial, PopulationVector{tt}, Para1, edges, Accu, PercentSheet);
    end
end

