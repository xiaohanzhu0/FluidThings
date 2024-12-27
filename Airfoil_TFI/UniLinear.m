clear
xy = geometry();
ineg = xy{1};
ipos = xy{2};
jneg = xy{3};
jpos = xy{4};

figure(1)
plot(ineg(:,1),ineg(:,2),'-o'); hold on
plot(ipos(:,1),ipos(:,2),'-o'); hold on
plot(jneg(:,1),jneg(:,2)); hold on
plot(jpos(:,1),jpos(:,2)); hold on
legend('ineg', 'ipos', 'jneg', 'jpos', Location='best');
%%
Nj = size(jneg, 1);
xi2 = linspace(0,1,Nj)';

Ni = 10;
xi1 = linspace(0,1,Ni)';

s = (1-xi1).*reshape(jneg,1,Nj,2) + xi1.*reshape(jpos,1,Nj,2);

%% Plot nodes
figure(2)
s1 = reshape(s(:,:,1),1,[]);
s2 = reshape(s(:,:,2),1,[]);
scatter(s1, s2, '.');
clear s1 s2

%% Save node points
save('UniLinear.mat','s');
%% Plot edges
figure(3)
for i = 1:Ni
    plot(s(i,:,1), s(i,:,2), 'r'); hold on
end

for j = 1:10:Nj % too dense, only keep 10% for plot
    plot(s(:,j,1), s(:,j,2), 'k'); hold on
end
plot(s(:,Nj,1), s(:,Nj,2), 'k'); hold on
