%一阶波动方程初边值问题的有限差分法
%齐次波动方程 齐次边界条件 初始条件为u(t=0)=sin(3pi*x/L),du/dt(t=0) = x(L-x) (0<x<L) 
%收敛条件 c=a*dt/dx <1
clear 
a = 2;
T =100; L = 4;  dt = 0.1;  dx = 0.3;
c = (a*dt/dx).^2;
x = linspace(0,L,ceil(L/dx));
%初始与边界条件  
u = zeros(length(x),ceil(T/dt));%每一列是某一时间所有位置的值，每一行是某一位置所有时间的值
u(:,1) = sin(3*pi*x/L);
u(:,2) = u(:,1)+(x.*(L-x)*dt)';
u(1,:) = 0;  u(end,:) = 0;%Dirichlet condition
h = plot(x,u(:,2),'linewidth',3);
axis ([0,L,-20,5]);
set(h,'EraseMode','xor')
for j = 2: ceil(T/dt)-1
    set(h,'XData',x,'YData',u(:,j));
    drawnow;
    pause(0.01)
    u(2:end-1,j+1) = c*u(3:end,j) + 2*(1-c)*u(2:end-1,j) +c*u(1:end-2,j)-u(2:end-1,j-1)-dt*0.1*(u(2:end-1,j)-u(2:end-1,j-1));%强迫振动项-dt*0.1*(u(2:end-1,j)-u(2:end-1,j-1))-dt.^2*9.8;
end
