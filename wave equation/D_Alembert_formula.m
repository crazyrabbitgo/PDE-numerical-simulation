%一维波动方程的柯西问题，达朗贝尔公式作图
%设初始位移为sin(7pai/L)x (3L/7<x<4L/7);此处设L=140
u = zeros(1,140);
x = linspace(0,1,140);
u(61:80) = 0.05*sin(pi*x(61:80)*7);
uu = u;
h = plot(x,u,'linewidth',3);
axis ([0,1,-0.05,0.05]);
set(h,'EraseMode','xor')

for at = 2:60
    lu(1:140) = 0;      ru(1:140) = 0 ;
    lx = [61:80] - at;  rx = [61:80]+at;
    lu(lx) = 0.5*uu(61:80); ru(rx) = 0.5*uu(61:80);
    u = lu + ru;
    set(h,'XData',x,'YData',u);
    drawnow;
    pause(0.2)
end