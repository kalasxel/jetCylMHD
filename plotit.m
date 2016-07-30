clear;
 
getDAT;

set(0,'DefaultAxesFontSize',10,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',18,'DefaultTextFontName','Times New Roman'); 

name = {'Rho', 'Vx', 'Vy', 'Vz', 'Bx', 'By', 'Bz', 'P'};
value = {'gm/cm^3', 'cm/c', 'cm/c', 'cm/c', 'Gs', 'Gs', 'Gs', 'dyn/cm^2'};

type = 1;

figure 
M=[];
for k = 1:3:size(TIME,2)-1
    
    txt = uicontrol('Style','text',...
        'Position',[400 45 120 20],...
        'String',[ 'time ' num2str(TIME(k),6) ]);
    
    subplot(2,2,1)
    type=1;
    plot(squeeze(GRID),squeeze(DATA(type,k,:)),'r','LineWidth',1.2)
    hold on
    title(name(type))
    xlabel(['radius, ' 'cm']);
    ylabel(value(type));
    grid on
    hold off
    
    subplot(2,2,2)
    type=2;
    plot(GRID,squeeze(DATA(type,k,:)),'r','LineWidth',1.2)
    hold on
    title(name(type))
    xlabel(['radius, ' 'cm']);
    ylabel(value(type));
    grid on
    hold off
    
    subplot(2,2,3)
    type=3;
    plot(GRID,squeeze(DATA(type,k,:)),'r','LineWidth',1.2)
    hold on
    title(name(type))
    xlabel(['radius, ' 'cm']);
    ylabel(value(type));
    grid on
    hold off
    
    subplot(2,2,4)
    type=4;
    plot(GRID,squeeze(DATA(type,k,:)),'r','LineWidth',1.2)
    hold on
    title(name(type))
    xlabel(['radius, ' 'cm']);
    ylabel(value(type));
    grid on
    hold off    
    M=[M,getframe(gcf)];
    
     pause(0.001)
    
end

%movie2avi(M,'test.avi');
% LeVeque R. J. Numer. Meth. for conv. laws