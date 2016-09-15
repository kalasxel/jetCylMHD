function PLOTall = PLOTall(type1,type2,type3,type4,gridX,gridY,RESULT,Time)


    nbLevels = 10;

    set(0,'DefaultAxesFontSize',10,'DefaultAxesFontName','Times New Roman');
    [x,y] = meshgrid(gridX,gridY);

    name = {'n', 'Vr', 'Vphi', 'Vz', 'Br', 'Bphi', 'Bz', 'T'};
    value = {'1/cm^3', 'cm/c', 'cm/c', 'cm/c', 'Gs', 'Gs', 'Gs', 'eV'};

    screenSize = get(0,'ScreenSize');     
    figure('ToolBar','none' , 'MenuBar','none')%, 'Position',screenSize);    
    M(1:size(Time,1)) = struct('cdata', [],'colormap', []);
    
    
    for k = 1:1:size(Time,1)
        colormap(hot)
        txt = uicontrol('Style','text',...
            'Position',[400 45 120 20],...
            'String',[ 'time ' num2str(Time(k),6) ]);

        subplot(2,2,1)
        type=type1;
            z = squeeze(RESULT(type,k,:,:));
            contourf(x,y,z,nbLevels,'LineColor','none');
        hold on
        colorbar
        title(name(type))
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off

        subplot(2,2,2)
        type=type2;
            z = squeeze(RESULT(type,k,:,:));
            contourf(x,y,z,nbLevels,'LineColor','none');
        hold on
        colorbar
        title(name(type))
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off

        subplot(2,2,3)
        type=type3;
            z = squeeze(RESULT(type,k,:,:));
            contourf(x,y,z,nbLevels,'LineColor','none');
        hold on
        colorbar
        title(name(type))
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off

        subplot(2,2,4)
        type=type4;
            z = squeeze(RESULT(type,k,:,:));
            contourf(x,y,z,nbLevels,'LineColor','none');
        hold on
        colorbar
        title(name(type))
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off 
        
        
        M(k)=getframe(gcf);

         pause(0.01)
        pause()

    end
% movie2avi(M,'test.avi');
end