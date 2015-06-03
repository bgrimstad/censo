
close all

A = [-4000 -3000;
       -60   -80];
b = [-100000 -2000]';
lb = [0 0]';
c = [-7000 -6000];
[x0,x1]=meshgrid(-1:0.1:26,-1:0.1:26);
f = -((7000 - 200.*x0).*x0 + (4000 - 140.*x1).*x1);
hold on
for i=-1:26,
    for j=-1:26,
        linOK = A*[i j]' >= b;
        posOK = i >= 0 && j >= 0;
        feasible = linOK(1) && linOK(2) && posOK;
        if feasible,
            plot(i,j,'w+');
        end
    end
end
plotregion(A,b,lb,[],[1 1 1],0)
set(gca,'Color',[0.7 0.7 0.7]);
contour(x0,x1,f,20);
xlabel('x_0 (Apples [tonnes])');
ylabel('x_1 (Bananas [tonnes])');
plot(16,12,'bo')
plot(0,0,'bo')
text(1,1,'z_0')
text(16.2,13,'x^*=(16,12,88640)')
legend(['Feasible points'])
axis ([-1 26 -1 26])
grid
