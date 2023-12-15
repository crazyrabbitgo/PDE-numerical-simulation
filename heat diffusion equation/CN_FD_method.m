function CN_FD_method
%�ó���Ϊʹ��CN��ʽ(��ʽ���)���һά�����ȴ������̵ĳ���ֵ����
clear 
global dx dt
dx = 0.1;%�ռ䲽��
dt = 0.1;%ʱ�䲽��
%��������
c = 3; p = 1; k = 1;%c,p,k�ֱ�Ϊ�����ݣ��ܶȺʹ���ϵ��
T_start = 0; T_end = 20; L_start = 0; L_end = 10; %T��LΪʱ��ͳ���
t = linspace(T_start,T_end,ceil((T_end-T_start)/dt));
x = linspace(L_start,L_end,ceil((L_end-L_start)/dx));

a = k/(c*p);%���̱�׼��ϵ��
lam = a*dt/(dx.^2);%��ַ��еĲ���

%���ƾ����ϵ������A
A = spdiags([-lam*ones(length(x),1) 1+2*lam*ones(length(x),1) ...
    -lam*ones(length(x),1)],[-1 0 1],length(x),length(x));
A_inv = inv(A);
u = zeros(length(t),length(x));%�洢���

L_fun = t.^0-1.1; R_fun = t.^0-1.1;%���ұ߽���������
%��ʼ����,���ȸ��ݳ�ʼ����ȷ��t=0ʱ��ֵ


u(1,:) = ic(x);
u(1,[1,end]) = [L_fun(1) R_fun(1)];%��t=0ʱ
f =@(x,t) x.^0+t.^0-2; %��Դ����,��ʱ�����0
for i = 2:length(t)
    ulur = bc(0,1,u(i-1,[1,end]),L_fun(i),R_fun(i));
    u(i,:) = Recur(A_inv,u(i-1,:),f(x,i*dt),ulur);
end
    function u_next = Recur(A_inv,u_now,fn,ulur)%�ú���Ϊ���ƺ�������������Ϊlength(x)
        u_next = A_inv*(u_now+fn.*dt)';
        u_next([1,end]) = ulur;%���ұ߽�����
    end

    function u0 = ic(x) %��ʼ������������Ϊ�ȴ�������ֻ��һ����ʼ����
        u0 = x.^0+9;
    end

    function ulur = bc(p,q,u_last,L_fun,R_fun)%�߽���������,��ʽΪp*u+q*dudx=fun 
        %u_last������һʱ������ұ߽�ֵ
        ul = q/(p*dx+q)*u_last(1)+dx/(p*dx+q)*L_fun;
        ur = q/(p*dx+q)*u_last(2)+dx/(p*dx+q)*R_fun;
        ulur = [ul ur];
    end

figure(1)%��������ͼ��
surf(x,t,u,'linestyle','none')
title('Thermal Conduction');
colormap('jet');
xlabel('Distance x');
ylabel('Time');

figure(2)
subplot(1,2,1)
[X,T] = meshgrid(x,t);%�����ȸ���ͼ��

colormap('hot')
contourf(X,T,u)
shading interp; 
colorbar
xlabel('Distance x');
ylabel('Time');

%�鿴�¶�����ͼ
n = 5;%��ʱ��T��Ϊn���׶Σ��鿴�����׶ν���ʱ���¶�����
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