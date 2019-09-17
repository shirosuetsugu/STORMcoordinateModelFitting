function simNoise = generateNoise3DtablePrint(FileName2, nPoints, box, plotflag)
% simNoise = generateNoise(nPoints, sizeVec)
% 
% generates Poisson noise with given number of points (nPoints) and given
% image size (sizeVec = [xlim1, xlim2, ylim1, ylim2, zlim1, zlim2]).
simNoise=table;

simNoise.x = box(1) + rand(nPoints,1)*(box(2) - box(1));
simNoise.y = box(3) + rand(nPoints,1)*(box(4) - box(3));
simNoise.z = box(5) + rand(nPoints,1)*(box(6) - box(5));

%% print
if plotflag==1
figure
subplot(2,2,1);
line([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)]);
hold on
scatter(simNoise.x,simNoise.y, '.');
axis equal

subplot(2,2,2);
line([box(1) box(1) box(2) box(2) box(1)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
scatter(simNoise.x,simNoise.z, '.');
axis equal

subplot(2,2,3);
line([box(3) box(3) box(4) box(4) box(3)], [box(5) box(6) box(6) box(5) box(5)]);
hold on
scatter(simNoise.y,simNoise.z), '.';
axis equal

subplot(2,2,4);
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(5) box(5) box(5) box(5) box(5)], 'Color', 'b');
hold on
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)], [box(6) box(6) box(6) box(6) box(6)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(3) box(3) box(3) box(3) box(3)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
plot3([box(1) box(1) box(2) box(2) box(1)], [box(4) box(4) box(4) box(4) box(4)], [box(5) box(6) box(6) box(5) box(5)], 'Color', 'b');
scatter3(simNoise.x,simNoise.y,simNoise.z, '.');
axis equal

[pathstr,name,ext] = fileparts(FileName2);

FigName=[pwd filesep name '-3D-NthNeighbor-Random-Simulated-'  num2str(nPoints) '-signals.pdf'];
print('-bestfit', FigName,'-dpdf','-r0');

end

end