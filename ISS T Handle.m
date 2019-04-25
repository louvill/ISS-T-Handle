%% Set up
clc
clear all

%% Problem 3 Part b
Ix = 0.00004942;
Iy = 0.00004317;
Iz = 0.00000901;
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];

tspan = [0 100];
icomega = [1; 0; 0; 0; pi/2; 0];

[t, omegab] = ode45(@omegaequations, tspan, icomega);

for i = 1:1:length(t)
    % H(i,:) = I*transpose(omegab(i,:));
    omegaatt = [omegab(i,1); omegab(i,2); omegab(i,3)];
    % if i == 3
    %     test = transpose(H(i-1,:))
    %     test2 = I*omegaatt
    % end
    H(:,i) = I*omegaatt;
    thetab(i) = acos(dot(H(:,i), omegaatt)/(norm(H(:,1))*norm(omegaatt)));
    % 4 is phi
    % 5 is theta
    % 6 is psi
    l = [cos(omegab(i,4))*cos(omegab(i,6))-cos(omegab(i,5))*sin(omegab(i,6))*sin(omegab(i,4))...
        cos(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*sin(omegab(i,6))*cos(omegab(i,4))...
        sin(omegab(i,5))*sin(omegab(i,6));
        -sin(omegab(i,6))*cos(omegab(i,4))-cos(omegab(i,5))*cos(omegab(i,6))*sin(omegab(i,4))...
        -sin(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*cos(omegab(i,4))*cos(omegab(i,6))...
        sin(omegab(i,5))*cos(omegab(i,6));
        sin(omegab(i,5))*sin(omegab(i,4)),-sin(omegab(i,5))*cos(omegab(i,4)),cos(omegab(i,5))];
    Hi(:,i) = transpose(l)*H(:,i);
    thetai(i) = acos(dot(Hi(:,i), omegaatt)/(norm(Hi(:,i))*norm(omegaatt)));
end

index = 1;
for i = 1:1:length(t)
    if imag(thetab(i)) == 0
        t2(index) = t(i);
        thetab2(index) = thetab(i);
    end
        
    index = index + 1;
end
figure(1)
hold on
%plot(t,transpose(thetab.*180./pi))
plot(transpose(t2),transpose(thetab2.*180./pi))
plot(t,transpose(thetai.*180./pi))
legend("Body System", "Inertial System")
grid on
xlabel("Time (s)",'Interpreter','Latex')
ylabel('Angle between H and $$\hat{\omega}$$ (deg)','Interpreter','Latex')
title('Angle between H and $$\hat{\omega}$$ with $$\omega_{x,t=0}$$ = 1 rad/s','Interpreter','Latex')

% figure(2)
% hold on
% plot(t, omegab(:,1))
% plot(t, omegab(:,2))
% plot(t, omegab(:,3))
% plot(t, omegab(:,4))
% plot(t, omegab(:,5))
% plot(t, omegab(:,6))
% legend("wx","wy","wz","phi","theta","psi")
% grid on
% xlabel("Time (s)")

%% Problem 3 Part c
tspan = [0 100];
icomega = [0; 1; 0.0001; 0; pi/2; 0];

[t, omegab] = ode45(@omegaequations, tspan, icomega);

for i = 1:1:length(t)
    % H(i,:) = I*transpose(omegab(i,:));
    omegaatt = [omegab(i,1); omegab(i,2); omegab(i,3)];
    % if i == 3
    %     test = transpose(H(i-1,:))
    %     test2 = I*omegaatt
    % end
    H(:,i) = I*omegaatt;
    thetab(i) = acos(dot(H(:,i), omegaatt)/(norm(H(:,1))*norm(omegaatt)));
    % 4 is phi
    % 5 is theta
    % 6 is psi
    l = [cos(omegab(i,4))*cos(omegab(i,6))-cos(omegab(i,5))*sin(omegab(i,6))*sin(omegab(i,4))...
        cos(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*sin(omegab(i,6))*cos(omegab(i,4))...
        sin(omegab(i,5))*sin(omegab(i,6));
        -sin(omegab(i,6))*cos(omegab(i,4))-cos(omegab(i,5))*cos(omegab(i,6))*sin(omegab(i,4))...
        -sin(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*cos(omegab(i,4))*cos(omegab(i,6))...
        sin(omegab(i,5))*cos(omegab(i,6));
        sin(omegab(i,5))*sin(omegab(i,4)),-sin(omegab(i,5))*cos(omegab(i,4)),cos(omegab(i,5))];
    Hi(:,i) = transpose(l)*H(:,i);
    thetai(i) = acos(dot(Hi(:,i), omegaatt)/(norm(Hi(:,i))*norm(omegaatt)));
end

index = 1;
for i = 1:1:length(t)
    if imag(thetab(i)) == 0
        if real(thetab(i)) == 0   
        else
            t2(index) = t(i);
            thetab2(index) = thetab(i); 
            index = index + 1;
        end
    end
end

figure(2)
hold on
%plot(t,transpose(thetab.*180./pi))
plot(transpose(t2),transpose(thetab2.*180./pi))
plot(t,transpose(thetai.*180./pi))
legend("Body System", "Inertial System")
grid on
xlabel("Time (s)",'Interpreter','Latex')
ylabel('Angle between H and $$\hat{\omega}$$ (deg)','Interpreter','Latex')
title('Angle between H and $$\hat{\omega}$$ with $$\omega_{y,t=0}$$ = 1 rad/s and $$\omega_{z,t=0}$$ = .0001 rad/s'...
    ,'Interpreter','Latex')
ylim([0 90])

%% Problem 3 Part d

% Video of part d: https://www.youtube.com/watch?v=0PL4awDoxtQ&feature=youtu.be
% Image files are attached
% Code commented for publishing of PDF

% % test = [1; 1; 1];
% % figure(3)
% hold on
% grid on
% % Thandle = [0 0 .05 0 -.05; .1/7+.02 -.04/7 -.04/7 -.04/7 -.04/7;...
% %     0 0 0 0 0];
% % plot3(0, .1/7+.02, 0, 'ko')
% % plot3(0, -.04/7, 0, 'ko')
% % plot3(.05, -.04/7, 0, 'ko')
% % plot3(-.05, -.04/7, 0, 'ko')
% % plot3(test(1,:), test(2,:), test(3,:))
% 
% 
% 
% for i = 1:1:length(t)
%     figure3 = figure(3)
%     hold off
%     point1 = [0; .1/7+.02; 0];
%     point2 = [0; -.04/7; 0];
%     point3 = [.05; -.04/7; 0];
%     point4 = [-.05; -.04/7; 0];
%     l = [cos(omegab(i,4))*cos(omegab(i,6))-cos(omegab(i,5))*sin(omegab(i,6))*sin(omegab(i,4))...
%         cos(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*sin(omegab(i,6))*cos(omegab(i,4))...
%         sin(omegab(i,5))*sin(omegab(i,6));
%         -sin(omegab(i,6))*cos(omegab(i,4))-cos(omegab(i,5))*cos(omegab(i,6))*sin(omegab(i,4))...
%         -sin(omegab(i,6))*sin(omegab(i,4))+cos(omegab(i,5))*cos(omegab(i,4))*cos(omegab(i,6))...
%         sin(omegab(i,5))*cos(omegab(i,6));
%         sin(omegab(i,5))*sin(omegab(i,4)),-sin(omegab(i,5))*cos(omegab(i,4)),cos(omegab(i,5))];
%     point1 = transpose(l)*point1;
%     point2 = transpose(l)*point2;
%     point3 = transpose(l)*point3;
%     point4 = transpose(l)*point4;
%     %plot3(point1(1),point1(2),point1(3),'ko')
%     %plot3([point1(1) point2(1)],[point1(2) point2(2)],[point1(3) point2(3)],'k-')
%     hold on
%     %plot3(point2(1),point2(2),point2(3),'ko')
%     %plot3(point3(1),point3(2),point3(3),'ko')
%     %plot3([point3(1) point4(1)],[point3(2) point4(2)],[point3(3) point4(3)],'k-')
%     %plot3(point4(1),point4(2),point4(3),'ko')
%     xlim([-.075 .075])
%     ylim([-.075 .075])
%     zlim([-.075 .075])
%     xlabel("X")
%     ylabel("Y")
%     zlabel("Z")
%     grid on
%     pause(.05)
%     %title = "THandle_" + i + ".png";
%     %saveas(figure3,title)
% end