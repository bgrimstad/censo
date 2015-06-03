
close all

A = [-4000 -3000;
       -60   -80];
b = [-100000 -2000]';
lb = [0 0]';
c = [-7000 -6000];
[x0,x1]=meshgrid(-1:0.1:26,-1:0.1:26);
f = -((7000 - 200.*x0).*x0 + (4000 - 140.*x1).*x1);
hold on
axis ([-1 26 -1 26])
for i=-1:26,
    for j=-1:26,
        if A*[i j]' >= b,
            plot(i,j,'w+');
        else
            plot(i,j,'k+');
        end
    end
end
plotregion(A,b,lb,[],[1 1 1],0)
set(gca,'Color',[0.8 0.8 0.8]);
contour(x0,x1,f,20);
xlabel('x_0 (Apples [tonnes])');
ylabel('x_1 (Bananas [tonnes])');
plot(15.7178,12.3762,'bo')
plot(0,0,'bo')
text(1,1,'z_0')
text(16.2,13,'x^*=(15.7178,12.3762,88676)')
legend('Feasible points')
grid
