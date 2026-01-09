clear
% parameters
Nx1 = 80;   Nx2 = 80;
s1  = linspace(0,1,Nx1);   ds1 = s1(2)-s1(1);
s2  = linspace(0,1,Nx2);   ds2 = s2(2)-s2(1);

w       = 1.7;    % SOR factor
w_ortho = 0.001;   % turn >0 to enable non-orth penalty
maxIter = 500;

% initial identity mesh
X = repmat(s1.',1,Nx2);
Y = repmat(s2, Nx1,1);

%[Y, X] = InitProb3(Nx1, Nx2, 0.1);
%boundary_points.b = [X(1,:); Y(1,:)];
%boundary_points.t = [X(end,:); Y(end,:)];
%boundary_points.l = [X(:,1), Y(:,1)];
%boundary_points.r = [X(:,end), Y(:,end)];



% --- INITIAL MESH (identity map) ---
X = repmat(s1.', 1, Nx2);   % X(i,j) = s1(i)
Y = repmat(s2,  Nx1, 1);    % Y(i,j) = s2(j)

for iter = 1:maxIter
  % (1) Update X along s1‐lines (constant j)
  %    compute metric at midpoints in s1
  Xmid = 0.5*(X(2:end,:)+X(1:end-1,:));    
  Ymid = 0.5*(Y(2:end,:)+Y(1:end-1,:));
  M11  = sqrt((1 + 15*Xmid).^(-2));
  %M11  = sqrt(1000 + 600*sin(2*pi*Xmid).*sin(2*pi*Ymid));

  for j = 1:Nx2
    for i = 2:(Nx1-1)
      x_old = X(i,j);
      % weighted average of neighbors in the s1‐direction
      t = ( M11(i-1,j)*X(i-1,j) + M11(i,j)*X(i+1,j) ) ...
          / ( M11(i-1,j) +   M11(i,j) );
      
      x1 = X(i-1,j); x2 = X(i+1,j);
      y1 = Y(i-1,j); y2 = Y(i+1,j); 
      if j==1; x3 = X(i,j); y3 = Y(i,j); x4 = X(i,j+1);  y4 = Y(i,j+1); 
      elseif j==Nx2;  x3 = X(i,j-1); y3 = Y(i,j-1); x4 = X(i,j); y4 = Y(i,j);
      else;  x3 = X(i,j-1); x4 = X(i,j+1);  y3 = Y(i,j-1); y4 = Y(i,j+1); end
      t_ortho = ((x1*y2 - y1*x2)*(x3-x4) - (x3*y4 - y3*x4)*(x1-x2)) / ((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4));
      t = t + w_ortho*t_ortho;


      X(i,j) = (1-w)*x_old + w*t;
    end
  end

  % (2) Update Y along s2‐lines (constant i)
  %    compute metric at midpoints in s2
  Xmid = 0.5*(X(:,2:end)+X(:,1:end-1));    
  Ymid = 0.5*(Y(:,2:end)+Y(:,1:end-1));           
  M22  = sqrt((1 + 15*Ymid).^(-2));  
  %M22 = sqrt(1000 - 600*sin(2*pi*Xmid).*sin(2*pi*Ymid));

  for i = 1:Nx1
    for j = 2:(Nx2-1)
      y_old = Y(i,j);
      % weighted average of neighbors in the s2‐direction
      t = ( M22(i,j-1)*Y(i,j-1) + M22(i,j)*Y(i,j+1) ) ...
          / ( M22(i,j-1) +   M22(i,j) );

      x3 = X(i,j-1); x4 = X(i,j+1); 
      y3 = Y(i,j-1); y4 = Y(i,j+1);
      if i==1; x1 = X(i,j); y1 = Y(i,j); x2 = X(i+1,j);  y2 = Y(i+1,j); 
      elseif i==Nx1;  x1 = X(i-1,j); y1 = Y(i-1,j); x2 = X(i,j);  y2 = Y(i,j); 
      else;  x1 = X(i-1,j); x2 = X(i+1,j);  y1 = Y(i-1,j); y2 = Y(i+1,j); end
      t_ortho = ((x1*y2 - y1*x2)*(y3-y4) - (x3*y4 - y3*x4)*(y1-y2)) / ((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4));
      t = t + w_ortho*t_ortho;

      Y(i,j) = (1-w)*y_old + w*t;
    end
  end

  % (optional) check convergence by max nodal move
  % maxMove = max( max(abs(X-Xprev)), max(abs(Y-Yprev)) );
  % if maxMove<tol, break, end

  % Xprev = X;  Yprev = Y;

  figure(1)
  plot(X,Y,'k'), hold on
  plot(X.',Y.','k'), hold off
  pause(0.2)


end

%% visualize final mesh
figure
plot(X,Y,'k.-'), hold on
plot(X.',Y.','k.-')
axis equal tight
xlabel('X'), ylabel('Y')
title('Converged 2D Grid via Alternating SOR')







%%
Nx1 = 40;
s1 = linspace(0,1,Nx1);
x1 = s1;
cost = [];

% --- FIX: Add an under-relaxation factor ---
w = 1.7;  % Start with a small value and increase if stable

for iter = 1:1000
    % The metric is calculated based on the state of x1 at the start of the iteration
    x1_int = (x1(2:end) + x1(1:end-1)) / 2;
    M_int = (1 + 100*15 * x1_int).^(-2);

    for i = 2:(Nx1-1)
        % Store the current position of the point before updating
        x_old_i = x1(i);

        % Calculate the target position based on neighbors and the (lagged) metric
        x_target = (sqrt(M_int(i-1))*x1(i-1) +  sqrt(M_int(i))*x1(i+1)) / (sqrt(M_int(i-1))+sqrt(M_int(i)));

        % --- FIX: Apply under-relaxation to the update ---
        % Move only part of the way to the target
        x1(i) = (1 - w) * x_old_i + w * x_target;
    end
    cost = [cost, sum(M_int.*(x1(2:end) - x1(1:end-1)).^2)];

    figure(1)
    plot(s1,x1)
    figure(2)
    plot((s1(2:end)+s1(1:end-1))/2,M_int.*(x1(2:end) - x1(1:end-1)).^2)
    figure(3)
    plot(cost)
    pause(0.1)
end





%%
Nx1 = 100;
s1 = linspace(0,1,Nx1);
x1 = s1;


for iter = 1:1000
    for i = 2:(Nx1-1)
        x1_int = (x1(2:end) + x1(1:end-1)) / 2;
        M_int = (1 + 1*15 * x1_int).^(-2);
        x1(i) = (M_int(i-1)*x1(i-1) +  M_int(i)*x1(i+1)) / (M_int(i-1)+M_int(i));
    end

    %for i = (Nx1-1):-1:2
    %    x1_int = (x1(2:end) + x1(1:end-1)) / 2;
    %    M_int = 40000 * (1 + 1*15 * x1_int).^(-2);
    %    x1(i) = (M_int(i-1)*x1(i-1) +  M_int(i)*x1(i+1)) / (M_int(i-1)+M_int(i));
    %end

    %x1_int = (x1(2:end) + x1(1:end-1)) / 2;
    %M_int = (1 + 1*15 * x1_int).^(-2);
    %even = 2:2:Nx1-1;
    %x1(even) = (M_int(even-1).*x1(even-1) +  M_int(even).*x1(even+1)) ./ (M_int(even-1)+M_int(even));

    %x1_int = (x1(2:end) + x1(1:end-1)) / 2;
    %M_int = (1 + 1*15 * x1_int).^(-2);
    %odd = 3:2:Nx1-1;
    %x1(odd) = (M_int(odd-1).*x1(odd-1) +  M_int(odd).*x1(odd+1)) ./ (M_int(odd-1)+M_int(odd));

    

    figure(1)
    plot(s1,x1)
    figure(2)
    plot((s1(2:end)+s1(1:end-1))/2,M_int.*(x1(2:end) - x1(1:end-1)))
    sum(M_int.*(x1(2:end) - x1(1:end-1)).^2)
    pause(0.1)
end