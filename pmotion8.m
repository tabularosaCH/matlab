function [] = pmotion8()
g = 9.81; m1=200; m2=10;x0=0; v0=0;
y0=0; yf=828; x0=0; xf=5; xdist = 5;
L = sqrt(xf^2+yf^2);

%PARAMETERS
A =5; %amplitude
w = sqrt(g/L); %angular speed
w_deg=sqrt(g/L)*(180/pi);
Tp = (2*pi)/w_deg; %period
freq = w_deg/(2*pi);

syms t X(t);
                    %disp('Dfferentiate X(t) wrt to t')
X = A*cos(w_deg*t); %X(t)
theta = xdist/L;
theta_deg =  theta*(180/pi);
u0=[theta_deg 0]; %initial angle, initial velocity

dxdt = 0.01;
tspan= 0:dxdt:1;

DX=@(t,u) [u(1); -((m1+m2)*g*sin(u(1)))/L]; %motion
[t,u] = ode45(DX, tspan,u0);

x = u(:,1); %displacement
v = u(:,2); %velocity
length(x), length(v)

%plot pendulum
figure(8)
clf; grid on; hold on; 
%%%%%%%%%%%%%%     
 
for i =1:length(t) 
    %plot pendulum  
    subplot(2,3,[1 2 4 5]) 
    
    xrect = [-2 2 2 1 0 -1 -2 -2]; %x values for building
    yrect = [-8.28 -8.28 -2 -2 0 -2 -2 -8.28];%y values for building
    plot(xrect,yrect) 
    set(gcf, 'Name', 'Simple Pendulum');
    axis(gca, 'equal'); %aspect ratio
    axis([-10 10 -12 2]);

    xlabel('x displacement'); ylabel('y displacement');
    title(sprintf(' Pendulum Animation\n  Time: %0.2f sec', t(i)));
    %%%%%%%%%
    text(-1,0.5,'Burj Khalifaa')
    text(0,-0.3,'\theta', 'FontSize',8,'FontName', 'Times')
    text(1.5,-4,'L')
    text(3,-8.3,sprintf('  x0\n(block)'))
    text(-3 ,-8.3,sprintf('  xp\n(block)'))
    %%%%%%%%%

    p1 = [0 0]; %position at equilibirum
    p2=  8.28*[sin(u(i,1)) -cos(u(i,1))]; %initial position at x0
    p3= 8.28*[sin(-theta_deg)  -cos(-theta_deg)];%peak position at xp

    circ_pend1 = viscircles(p2, 0.2); %block at x0
    circ_pend2 = viscircles(p3, 0.2, 'LineStyle',':'); %block at xp

    pendbar1 = line([p1(1) p2(1)], [p1(2) p2(2)]); %string at x0 
    pendbar2 = line([p1(1) p3(1)], [p1(2) p3(2)],'LineStyle','--');%string at xp 
    drawnow
    pause(0.001); %

    %plot displacement
    subplot(2,3,[3]);
    plot(t,x)
    xlabel('time [s]'), ylabel('displacement [m]'),title('Displacement vs. Time')
    
    %plot velocity
    subplot(2,3,[6]);
    plot(t,v)
    xlabel('time [s]'), ylabel('velocity [m/s]'),title('Velocity vs. Time')
    hold off;
end%for

end%pmotion8

