function one_dim_heat
%һά�ȴ�������,ʹ��pdepe���
x = linspace(0,pi,40);
t = linspace(0,10,30);
m = 0;

sol = pdepe(m,@heatpde,@heatic,@heatbc,x,t);

surf(x,t,sol(:,:,1),'linestyle','none')
title('Thermal Conduction');
colormap('hot');
xlabel('Distance x');
ylabel('Time');

[X,T] = meshgrid(x,t);
figure()
colormap('hot')
contourf(X,T,sol(:,:,1))
shading interp; 
colorbar; colormap(jet);

result = sol(:,:,1);


    function [c,f,s] = heatpde(x,t,u,dudx) % Equation to solve
        c = 1;
        f = dudx;
        s = cos(x/2);
    end
%----------------------------------------------
    function u0 = heatic(x) % Initial conditions
        u0 = x*cos(x/2);
    end
%----------------------------------------------
    function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t) % Boundary conditions
        pl = -1;%��߽纯����
        ql = 1;%��߽絼����ϵ��
        pr = ur-pi;%�ұ߽纯����
        qr = 0;%�ұ߽絼����ϵ��
    end
end