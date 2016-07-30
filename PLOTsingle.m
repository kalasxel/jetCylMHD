function PLOTsingle = PLOTsingle(type,gridX,gridY,RESULT,Time)
    set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
    [x,y] = meshgrid(gridX,gridY);

    name = {'Rho', 'Vr', 'Vphi', 'Vz', 'Br', 'Bphi', 'Bz', 'Ptot', 'P', 'e'};
    value = {'gm/cm^3', 'cm/c', 'cm/c', 'cm/c', 'Gs', 'Gs', 'Gs', 'dyn/cm^2', 'dyn/cm^2', 'dyn/cm^2'};

    screenSize = get(0,'ScreenSize'); 
    figure('ToolBar','none' , 'MenuBar','none', 'Position',screenSize);
%     writerObj = VideoWriter('peaks.avi');
%     open(writerObj);
    
    for k = 1:1:size(Time,1)
    colormap(hot)
        txt = uicontrol('Style','text',...
           'FontName', 'Arial', 'FontSize', 10, ...
          'Position',[400 45 120 20],...
          'String',[ 'time ' num2str(Time(k),6) ]);
    
        z = squeeze(RESULT(type,k,:,:));
        contourf(x,y,z,20,'LineColor','none');
        colorbar('vert')
        hold on
%         colorbar('Limits',[0,2],'LimitsMode','manual');
        title(name(type))
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off
        
%         frame = getframe;
%    writeVideo(writerObj,frame);
    
%         pause(0.00001)
          pause()
    
    end

%     movie2avi(M,'test.avi');
close(writerObj);
end
