clear
BCtype = 1; % 1 for linear, 2 for curved
BCspace = 2; % 1 for equi-spaced, 2 for non equi-spaced
xy = geometry();
ineg = xy{1};
ipos = xy{2};
jneg = xy{3};
jpos = xy{4};

%% Make equispaced data points in i boundary
Ni = 10;
Nj = size(jneg, 1);
[ineg, ipos] = make_boundary_data(jneg, jpos, Ni, BCspace);

if BCtype == 1
    ineg(2:end-1,1) = ineg(2:end-1,1);
elseif BCtype == 2
    ineg(2:end-1,1) = ineg(2:end-1,1) - 0.3*sin(ineg(2:end-1,2)*pi/1.5);
end

figure(1)
plot(ineg(:,1),ineg(:,2),'-o'); hold on
plot(ipos(:,1),ipos(:,2),'-o'); hold on
plot(jneg(:,1),jneg(:,2)); hold on
plot(jpos(:,1),jpos(:,2)); hold on
legend('ineg', 'ipos', 'jneg', 'jpos', Location='best');

%%
P12 = ineg(1,:);
P34 = ipos(end,:);
P14 = ineg(end,:);
P32 = ipos(1,:);

xi1 = linspace(0,1,Ni)';
xi2 = linspace(0,1,Nj)';

[XI1, XI2] = meshgrid(xi1, xi2);
s1 = (1-xi2)*ineg(:,1)' + xi2*ipos(:,1)' + ((1-xi1)*jneg(:,1)' + xi1*jpos(:,1)')' - ...
    ( (1-xi1)*(1-xi2')*P12(1) + xi1*xi2'*P34(1) + xi1.*(1-xi2')*P14(1) + (1-xi1).*xi2'*P32(1) )';
s1 = s1';

s2 = (1-xi2)*ineg(:,2)' + xi2*ipos(:,2)' + ((1-xi1)*jneg(:,2)' + xi1*jpos(:,2)')' - ...
    ( (1-xi1)*(1-xi2')*P12(2) + xi1*xi2'*P34(2) + xi1.*(1-xi2')*P14(2) + (1-xi1).*xi2'*P32(2) )';
s2 = s2';

%% Plot edges
figure(3)
for i = 1:Ni
    plot(s1(i,:), s2(i,:), 'r'); hold on
end

for j = 1:10:Nj % too dense, only keep 10% for plot
    plot(s1(:,j), s2(:,j), 'k'); hold on
end
plot(s1(:,Nj), s2(:,Nj), 'k'); hold on

%% Save node points
%save('BiLinear.mat','s1','s2');

function [ineg, ipos] = make_boundary_data(jneg, jpos, Ni, BCspace)
    if BCspace == 1
        Ix_neg = linspace(jneg(1,1), jpos(1,1), Ni);
        Iy_neg = linspace(jneg(1,2), jpos(1,2), Ni);
        Ix_pos = linspace(jneg(end,1), jpos(end,1), Ni);
        Iy_pos = linspace(jneg(end,2), jpos(end,2), Ni);
    elseif BCspace == 2
        cheb = -cos(linspace(0,pi,Ni)); % Use Chebyshev points for example
        Ix_neg = (jneg(1,1)+ jpos(1,1))/2 + (jpos(1,1)-jneg(1,1))*cheb/2;
        Iy_neg = (jneg(1,2)+ jpos(1,2))/2 + (jpos(1,2)-jneg(1,2))*cheb/2;
        Ix_pos = (jneg(end,1)+ jpos(end,1))/2 + (jpos(end,1)-jneg(end,1))*cheb/2;
        Iy_pos = (jneg(end,2)+ jpos(end,2))/2 + (jpos(end,2)-jneg(end,2))*cheb/2;
    end
    
    ineg = [Ix_neg', Iy_neg'];
    ipos = [Ix_pos', Iy_pos'];
end