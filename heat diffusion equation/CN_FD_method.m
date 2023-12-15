function CN_FD_method
%该程序为使用CN格式(隐式差分)求解一维线性热传导方程的初边值问题
clear 
global dx dt
dx = 0.1;%空间步长
dt = 0.1;%时间步长
%参数设置
c = 3; p = 1; k = 1;%c,p,k分别为比热容，密度和传热系数
T_start = 0; T_end = 20; L_start = 0; L_end = 10; %T和L为时间和长度
t = linspace(T_start,T_end,ceil((T_end-T_start)/dt));
x = linspace(L_start,L_end,ceil((L_end-L_start)/dx));

a = k/(c*p);%方程标准型系数
lam = a*dt/(dx.^2);%差分法中的参数

%递推矩阵的系数矩阵A
A = spdiags([-lam*ones(length(x),1) 1+2*lam*ones(length(x),1) ...
    -lam*ones(length(x),1)],[-1 0 1],length(x),length(x));
A_inv = inv(A);
u = zeros(length(t),length(x));%存储结果

L_fun = t.^0-1.1; R_fun = t.^0-1.1;%左右边界条件函数
%开始计算,首先根据初始条件确定t=0时的值


u(1,:) = ic(x);
u(1,[1,end]) = [L_fun(1) R_fun(1)];%当t=0时
f =@(x,t) x.^0+t.^0-2; %热源函数,此时恒等于0
for i = 2:length(t)
    ulur = bc(0,1,u(i-1,[1,end]),L_fun(i),R_fun(i));
    u(i,:) = Recur(A_inv,u(i-1,:),f(x,i*dt),ulur);
end
    function u_next = Recur(A_inv,u_now,fn,ulur)%该函数为递推函数，向量长度为length(x)
        u_next = A_inv*(u_now+fn.*dt)';
        u_next([1,end]) = ulur;%左右边界修正
    end

    function u0 = ic(x) %初始条件函数，因为热传导方程只有一个初始条件
        u0 = x.^0+9;
    end

    function ulur = bc(p,q,u_last,L_fun,R_fun)%边界条件函数,格式为p*u+q*dudx=fun 
        %u_last代表上一时间的左右边界值
        ul = q/(p*dx+q)*u_last(1)+dx/(p*dx+q)*L_fun;
        ur = q/(p*dx+q)*u_last(2)+dx/(p*dx+q)*R_fun;
        ulur = [ul ur];
    end

figure(1)%画出函数图像
surf(x,t,u,'linestyle','none')
title('Thermal Conduction');
colormap('jet');
xlabel('Distance x');
ylabel('Time');

figure(2)
subplot(1,2,1)
[X,T] = meshgrid(x,t);%画出等高线图像

colormap('hot')
contourf(X,T,u)
shading interp; 
colorbar
xlabel('Distance x');
ylabel('Time');

%查看温度曲线图
n = 5;%将时间T分为n个阶段，查看各个阶段结束时的温度曲线
subplot(1,2,2)
plot(x,u(1,:))
hold on
plot(x,u(ceil([1:T_end/n:T_end]/dt),:))

leg_str{1} = ['T =',num2str(0)];  
for i = 1:n
    leg_str{i+1} = ['T =',num2str(T_end/5*i)];  
end
legend(leg_str)
end