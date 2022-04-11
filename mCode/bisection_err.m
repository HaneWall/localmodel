function [z] = bisection_err(f, x_min, x_max, err)
    x_mid = 1/2 * (x_min + x_max);
    while abs(f(x_mid))>=err
        if f(x_mid)<0 && f(x_min)<0
            x_min = x_mid;
        else
            x_max = x_mid;
        end
        x_mid = 1/2 * (x_min + x_max);
    end
    z = x_mid;
end

