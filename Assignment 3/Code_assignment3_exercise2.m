function Chouliaras_assignment3_exercise2
% G.C. Chouliaras
%This function determines the center of mass of a domain enclosed by an
%interpolated curve.
%There is no input or output.
%The function asks the user to put his hand on a rectangle on the screen
%and outline the circumference of his hand using mouse clicks. The user can
%either press 'Yes' and proceed or press 'Cancel' to cancel the function.
%The function produces a figure of the input points, the interpolated curve
%and the area inside the curve, along with a mark of the center of mass.
%
%
%Example syntax:
% Chouliaras_assignment3_exercise2

button = questdlg('Put your hand on the rectangle on the screen and outline the circumference of the hand using mouse clicks.',...
    'Instructions','Yes','Cancel','Yes');

switch button
    case 'Yes'
        figure (1);
        FigHandle = figure('Position', [100, 100, 1024, 896]);
        hold on;
        [x,y] = ginput;
        hold off;
        %check for identical subsequent points
        counter = 1;
        while counter < length(x);
            if (x(counter) == x(counter+1)) && (y(counter) == y(counter+1))
                x(counter) = [];
                y(counter) = [];
            else
                counter = counter + 1;
            end
            
        end

        close Figure 1;

        %arclength parametrization
        
        
        n = length(x);
        %close the curve
        x(n+1) = x(1);
        y(n+1) = y(1);
        
        t = zeros(length(x),1);
        for index = 1:(length(x)-1)
            t(index+1) = t(index) + sqrt((x(index+1)-x(index))^2+(y(index+1) - y(index))^2);
        end
           
        %create a vector of 800 points between t(1) and t(n)
        tvector = linspace(t(1), t(n),800);
        
        
        %fit the splines
        
        %put spline throught (t,x)
        xspline = ppval(spline(t,x),tvector);
        
        %put spline through (t,y)
        
        yspline = ppval(spline(t,y),tvector);
        
        
        
        %computes the integrals using trapz
        
        
        integral_x = abs(trapz(yspline, xspline.^2/2));
        integral_y = abs(trapz(xspline, -yspline.^2/2));
        area = abs(trapz(yspline,xspline));
        
        %use the integrals to get the coordinates of the center of mass
        x_coordinate = integral_x /area;
        y_coordinate = integral_y/area;
        
        %plot the center of mass, the input points, the area under the
        %curve and the interpolated curve
        figure(2);
        plot(xspline,yspline,'k')
        fill(xspline,yspline,'c');
        hold on;
        plot(x_coordinate,y_coordinate,'kx','MarkerSize',15)
        hold on;
        plot(x,y,'r*')
        axis([0 1 0 1])
        hold off;
        
        
    case 'Cancel'
        return;
end


