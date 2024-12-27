function [xy] = geometry()
%
%  Define the boundaries
%  Must be defined in increasing index direction
%

%  airfoil and farfield
load airfoil_small
xy_jneg = flipud(xy_airfoil);
xy_jpos = flipud(xy_farfield);
clear xy_airfoil xy_farfield

%  close the airfoil
xy_jneg(1,:) = xy_jneg(end,:);

%  outflow
xy_ineg = [xy_jneg(1,:); xy_jpos(1,:)];
xy_ipos = [xy_jneg(end,:); xy_jpos(end,:)];

xy = {xy_ineg, xy_ipos, xy_jneg, xy_jpos};

%  check corners
err(1) = max(abs(xy_ineg(1,:) - xy_jneg(1,:)));
err(2) = max(abs(xy_ineg(end,:) - xy_jpos(1,:)));
err(3) = max(abs(xy_ipos(1,:) - xy_jneg(end,:)));
err(4) = max(abs(xy_ipos(end,:) - xy_jpos(end,:)));
if max(err) > 1e-8,
  disp('Error in geometry definition')
  pause
end


return
figure(10);clf
col='bcmr';
for c=1:4,
  plot(xy{c}(:,1), xy{c}(:,2), col(c), 'linewidth',2); hold on
end
legend('ineg','ipos','jneg','jpos')






