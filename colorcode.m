function c = colorcode(N)
if mod(N,7)==0
    c=[1 0 0];%'red'
elseif mod(N,7)==1
    c=[1 1 0];%'yellow'
elseif mod(N,7)==2
    c=[1 0 1];%'magenta'
elseif mod(N,7)==3
    c=[0 1 0];%'green'
elseif mod(N,7)==4
    c=[0 0 1];%'blue'
elseif mod(N,7)==5
    c=[0 1 1];%'cyan'
elseif mod(N,7)==6
    c=[0 0 0];%'black'

end
