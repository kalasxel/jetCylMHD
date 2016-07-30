clear

x_ini=0.0;
x_fin=1.0;
dx=0.002;
t_fin=1;
dt=0.0001;

nbR=int16((x_fin-x_ini)/dx) +1;
nbT=int16(t_fin/dt) +1;
clear x_ini x_fin dx t_fin dt

strR = '';
for i=1:nbR
    %strR = [strR  ' ' num2str(i)];
    strR = [strR  ' %f'];
end

RERO='single_results/0_Ro.dat';
REVR='single_results/1_Vx.dat';
REVP='single_results/2_Vy.dat';
REVZ='single_results/3_Vz.dat';
REBR='single_results/4_Bx.dat';
REBP='single_results/5_By.dat';
REBZ='single_results/6_Bz.dat';
% REPTOT='results/8_Ptot.dat';

GRID='single_results/GRID.dat';
TIME='single_results/TIME.dat';


DATA(1) = fopen(RERO,'r');
DATA(2) = fopen(REVR,'r');
DATA(3) = fopen(REVP,'r');
DATA(4) = fopen(REVZ,'r');
DATA(6) = fopen(REBP,'r');
DATA(7) = fopen(REBZ,'r');
% DATA(8) = fopen(REPTOT,'r');

SPATIAL = fopen(GRID,'r');
LAYER = fopen(TIME,'r');

clear RERO REVR REVP REVZ REBP REBZ REPTOT GRID TIME


RO=textscan(DATA(1),strR,nbT+5);
VR=textscan(DATA(2),strR,nbT+5);
VP=textscan(DATA(3),strR,nbT+5);
VZ=textscan(DATA(4),strR,nbT+5);
BP=textscan(DATA(6),strR,nbT+5);
BZ=textscan(DATA(7),strR,nbT+5);
% PT=textscan(DATA(8),strR,nbT+5);

SP=textscan(SPATIAL,'%f',nbR+1);
LA=textscan(LAYER,'%f',nbT+1);

err_cl=fclose('all');

clear DATA SPATIAL LAYER strR err_cl


for r=1:nbR
    for t=1:size(RO{1},1)
        DATA(1,t,r) = RO{r}(t);
        DATA(2,t,r) = VR{r}(t);
        DATA(3,t,r) = VP{r}(t);
        DATA(4,t,r) = VZ{r}(t);
        DATA(6,t,r) = BP{r}(t);
        DATA(7,t,r) = BZ{r}(t);
%         DATA(8,t,r) = PT{r}(t);    
        %TEST(t,r) = RO{r}(t);
    end  
end

for r=1:size(SP{1},1)
    GRID(r)=SP{1}(r);
end

for t=1:size(LA{1},1)
    TIME(t)=LA{1}(t);
end;

clear RO VR VP VZ BP BZ PT SP LA nbR nbT
clear i r t