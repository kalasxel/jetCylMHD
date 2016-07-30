clear

nbF = 10;

x_ini=0.0;
% x_fin=6.283185307179586;
% dx=0.049087385212340;
x_fin=5;
dx=0.05;    
NX=int16((x_fin-x_ini)/dx)+1;
y_ini=0.0;
% y_fin=6.283185307179586;
% dy=0.049087385212340;
y_fin=40;
dy=0.05;
NY=int16((y_fin-y_ini)/dy)+1;
t_fin=1;
dt=0.0001;
NT = int16(t_fin/dt)+1;

clear x_ini x_fin dx y_ini y_fin dy t_fin dt




type=cell(nbF,1);
type{1} = 'Ro';
type{2} = 'Vx';
type{3} = 'Vy';
type{4} = 'Vz';
type{5} = 'Bx';
type{6} = 'By';
type{7} = 'Bz';
type{8} = 'E';
type{9} = 'P';
type{10} = 'Ptot';
res='results/';

directory=cell(nbF,1);
for f=1:nbF
    directory{f} = cellstr(strcat(res, num2str(f-1), '_', type(f), '/'));
    %     disp(directory{f});
end

number = size(dir(char(directory{1})),1)-2;

names=cell(nbF,number);

for f=1:nbF
    for i=1:number
        names{f,i} = strcat(directory{f}, type{f}, '_', num2str(i-1), '.dat');
    end
end

strX='';
for i=1:NX
    strX = [strX  ' %f'];
end

RESULT = zeros(nbF,number,NY,NX);

for f=1:nbF
    for i=1:number
        
        IDopen = fopen(char(names{f,i}), 'r');
        
        if IDopen==-1
            disp(strcat('CANNOT OPEN FILE ', char(names{f,i}) ));
        end 
        
        readFile = textscan(IDopen, strX, NT+1);
        
        closeErr=fclose(IDopen);
        if closeErr==-1
            disp(strcat('CANNOT CLOSE FILE ', char(names{f,i}) ));
        end
        
        for n=1:NX
            for m=1:NY
                RESULT(f,i,m,n) = readFile{n}(m);
            end
        end
        
    end
%     disp(strcat(type{f}, ' has been read'));
end

clear directory closeErr i IDopen m n names readFile strX type


gridType = cell(3,1);
gridType{1} = strcat(res, 'GRIDR.dat');
gridType{2} = strcat(res, 'GRIDZ.dat');
gridType{3} = strcat(res, 'TIME.dat');

for f=1:3
    
    IDopen = fopen(char(gridType{f}), 'r');
    if IDopen==-1
        disp(strcat('CANNOT OPEN FILE ', char(gridType{f}) ));
    end    
    readFile(f) = textscan(IDopen, '%f');
    closeErr=fclose(IDopen);
    if closeErr==-1
        disp(strcat('CANNOT CLOSE FILE ', char(gridType{f}) ));
    end   
    
    for i=1:size(readFile{f},1)
        AllGrid(i,f) = readFile{f}(i);
    end    
    
end

gridX=AllGrid(1:NX,1);
gridY=AllGrid(1:NY,2);
Time=AllGrid(1:size(readFile{3},1),3);

clear AllGrid closeErr f i IDopen number readFile gridType res NX NY NT nbF

% save ModelResults.mat RESULT Time gridX gridY 
