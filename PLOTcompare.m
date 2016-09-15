function PLOTall = PLOTall(type1,type2,gridX,gridY,RESULT1,RESULT2,Time)


    set(0,'DefaultAxesFontSize',10,'DefaultAxesFontName','Times New Roman');
    [x,y] = meshgrid(gridX,gridY);

    name = {'n', 'Vr', 'Vphi', 'Vz', 'Br', 'Bphi', 'Bz', 'T'};
    value = {'1/cm^3', 'cm/c', 'cm/c', 'cm/c', 'Gs', 'Gs', 'Gs', 'eV'};

    screenSize = get(0,'ScreenSize');     
    figure('ToolBar','none' , 'MenuBar','none')%, 'Position',screenSize./2);    
    M(1:300) = struct('cdata', [],'colormap', []);

% writerObj = VideoWriter('peaks.avi');
% open(writerObj);

    for k = 1:1:300
        colormap(hot)
        
        txt = uicontrol('Style','text',...
            'Position',[400 45 120 20],...
            'String',[ 'time ' num2str(Time(k),6) ]);

        subplot(1,2,1)
        type=type1;
            z = squeeze(RESULT1(type,k,:,:));
            contourf(x,y,z,10,'LineColor','none');

        hold on
        colorbar
        title([name(type) 'RESULT1' ])
        xlabel(['radius, ' 'cm']);
        ylabel(value(type));
        grid on
        hold off

        subplot(1,2,2)
        type=type2;
            z = squeeze(RESULT2(type,k,:,:));
            contourf(x,y,z,10,'LineColor','none');
        hold on
        colorbar
        title([name(type) 'RESULT2'])
        xlabel(['radius, ' 'cm']);
        ylabel([value(type)]);
        grid on
        hold off

        M(k)=getframe(gcf);
%    frame = getframe;
%    writeVideo(writerObj,frame);


         pause(0.04)
%         pause()

    end
movie2avi(M,'test.avi');
% close(writerObj);
end