function [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2)
    dx2ds2 = zeros(size(x2));
    dx2ds2(2:end-1,:) = (x2(3:end,:)-x2(1:end-2,:)) / (sigma2*2);
    dx2ds2(1,:) = (-x2(3,:)+4*x2(2,:)-3*x2(1,:)) / (sigma2*2);
    dx2ds2(end,:) = -(-x2(end-2,:)+4*x2(end-1,:)-3*x2(end,:)) / (sigma2*2);
    %dx2ds2(1,:) = (x2(2,:)-x2(1,:)) / (sigma2);
    %dx2ds2(end,:) = (x2(end,:)-x2(end-1,:)) / (sigma2);

    dx1ds1 = zeros(size(x1));
    dx1ds1(:,2:end-1) = (x1(:,3:end)-x1(:,1:end-2)) / (sigma1*2);
    dx1ds1(:,1) = (-x1(:,3)+4*x1(:,2)-3*x1(:,1)) / (sigma1*2);
    dx1ds1(:,end) = -(-x1(:,end-2)+4*x1(:,end-1)-3*x1(:,end)) / (sigma1*2);
    %dx1ds1(:,1) = (x1(:,2)-x1(:,1)) / (sigma2);
    %dx1ds1(:,end) = (x1(:,end)-x1(:,end-1)) / (sigma2);
end