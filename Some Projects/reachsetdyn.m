function reachsetdyn(t1, t2, N, filename)
    showTrace = 1;
    showCaption = 0;
    whole = 0;
    showLegend = 1;
    
    if (showTrace == 1)
        showLegend = 0;
    end
    
    xMin = -2;
    xMax = 2;
    yMin = -6;
    yMax = 3;
    
%    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig = figure();
    l1 = plot(0, 0);
    l1.Color = 'red';
    hold on;
    l2 = plot(0, 0);
    l2.Color = 'black';
    if (showLegend == 1)
     legend('Граница множества достижимости', ...
            'Кривая первого переклюения');
    end
    h = (t2 - t1) / N;
    axis([xMin xMax yMin yMax]);
    xlabel('x_1');
    ylabel('x_2');
    if (filename ~= ' ')
        %% Initialize video
        fig.Visible = 'off';
        filename = strcat(filename, '.avi');
        myVideo = VideoWriter(filename); %open video file
        myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
        open(myVideo)
        %% Plot in a loop and grab frames
        for i = 1 : N
            [X, Y, Sx, Sy] = reachset(t1 + h * i, whole, showTrace);
            l1.XData = X;
            l1.YData = Y;
            l2.XData = Sx;
            l2.YData = Sy;
            if (showCaption == 1)
                t = strcat('Граница множества достижимости, T = ', ...
                    num2str(t1 + h * i));
                title(t);
            end
            pause(0.05); %Pause and grab frame
            frame = getframe(gcf); %get frame
            writeVideo(myVideo, frame);
        end
        close(myVideo);
        disp('Запись закончена');
    else
        for i = 1 : N
            [X, Y, Sx, Sy] = reachset(t1 + h * i, whole, showTrace);
            pause(0.05);
            l1.XData = X;
            l1.YData = Y;
            l2.XData = Sx;
            l2.YData = Sy;
            if (showCaption == 1)
                t = strcat('Граница множества достижимости, T = ', ...
                    num2str(t1 + h * i));
                title(t);
            end
        end
    end
end

















