%һ�ײ������̳���ֵ��������޲�ַ�
%��β������� ��α߽����� ��ʼ����Ϊu(t=0)=sin(3pi*x/L),du/dt(t=0) = x(L-x) (0<x<L) 
%�������� c=a*dt/dx <1
clear 
a = 2;
T =100; L = 4;  dt = 0.1;  dx = 0.3;
c = (a*dt/dx).^2;
x = linspace(0,L,ceil(L/dx));
%��ʼ��߽�����  
u = zeros(length(x),ceil(T/dt));%ÿһ����ĳһʱ������λ�õ�ֵ��ÿһ����ĳһλ������ʱ���ֵ
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
    u(2:end-1,j+1) = c*u(3:end,j) + 2*(1-c)*u(2:end-1,j) +c*u(1:end-2,j)-u(2:end-1,j-1)-dt*0.1*(u(2:end-1,j)-u(2:end-1,j-1));%ǿ������-dt*0.1*(u(2:end-1,j)-u(2:end-1,j-1))-dt.^2*9.8;
end
