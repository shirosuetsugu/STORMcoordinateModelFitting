function RandomDist3D=generateRandom3(NumSimData,Dist, Accu)
RandomDist3D=table;
RandomDist3D.x=zeros(NumSimData,1);
RandomDist3D.y=zeros(NumSimData,1);
RandomDist3D.z=zeros(NumSimData,1);

% add the additional data with gaussian dist. Dist variance = Accu
RandomDist3D.y=RandomDist3D.y+Dist+Accu*randn([NumSimData 1]);

for m=1:NumSimData
    % rotate dimers at Z axis
    tZ=(2*rand-1)*pi();
    Rz = [cos(tZ) -sin(tZ) 0; sin(tZ) cos(tZ) 0; 0 0 1];
    % rotate ximers at X axis
    tX=(2*rand-1)*pi();
    Rx = [1 0 0; 0 cos(tX) -sin(tX); 0 sin(tX) cos(tX)];
    xyz=[RandomDist3D.x(m);RandomDist3D.y(m);RandomDist3D.z(m)];
    xyzRx = Rz * (Rx * xyz);
    
    RandomDist3D.x(m)=xyzRx(1);
    RandomDist3D.y(m)=xyzRx(2);
    RandomDist3D.z(m)=xyzRx(3);
    
end
end
