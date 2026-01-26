function metric_data = metric_transform_inverse(metric_data, params)
%METRIC_TRANSFORM_INVERSE Compute inverse metric for Mp.
    detJ = metric_data.Mp11.*metric_data.Mp22 - metric_data.Mp12.^2;
    metric_data.Mp_inv11 = metric_data.Mp22 ./ detJ;
    metric_data.Mp_inv22 = metric_data.Mp11 ./ detJ;
    metric_data.Mp_inv12 = -metric_data.Mp12 ./ detJ;

    if params.diagonal_reg == 1
        metric_data.Mp_inv11 = metric_data.Mp_inv11 + 0.01*max(metric_data.Mp_inv11,[],'all');
        metric_data.Mp_inv22 = metric_data.Mp_inv22 + 0.01*max(metric_data.Mp_inv22,[],'all');
    end
end
