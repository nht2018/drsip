%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-12 20:56:52
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-16 21:27:39
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function x = prox_box(z, mu, l, u, implement)
    % compute the proximal operator of the log barrier of box constraint
    % i.e., solve x \in [l, u] from the equation
    % x - z - \mu / ( x - l) + \mu / (u - x) = 0
    % where z is the input, mu is the barrier parameter, l and u are the lower and upper bounds

    if  nargin < 5
        implement = 2;
    end
    assert(length(z) == length(l) && length(z) == length(u));
    assert(all(l <= u));

    if implement == 1
        %% implement 1: use root function. this is slower since we need a loop. 
        x = zeros(size(z));
        % Coefficients of the polynomial
        a = -1;
        b = l + u + z;
        c = - u.* z - l.* z - l.* u + 2 * mu;
        d = u .* l .* z - mu .* (l + u);
        
        for i = 1:length(z)
            % Define the polynomial
            p = [a b(i) c(i) d(i)];

            % Find the roots
            roots_p = roots(p);
            roots_p = real(roots_p); 

            % Filter roots that are within the interval [l(i), u(i)]
            x(i) = roots_p(roots_p >= l(i) & roots_p <= u(i));

        end
    elseif implement == 2
        % implement 2: use the root formula for cubic equation. this is faster but complex number may appear
        % Calculate the components of the solution
        % the following expression is from mathematica
        A = -((l + u - 2*z) .* (2*l.^2 - 9*mu + 2*u.^2 + u.*z - z.^2 + l.*(-5*u + z)));
        B = -(l + u + z).^2 - 3*(2*mu - l.*u - l.*z - u.*z);
        C = sqrt(A.^2 + 4*B.^3);

        % Constants for the cubic roots of unity
        omega1 = exp(2i * pi / 3);
        omega2 = exp(4i * pi / 3);
        
        temp = (A + C).^(1/3);
        x = (l + u + z)/3 + omega1 * (2^(1/3) * B) ./ (3 * temp) - omega2 * temp ./ (3 * 2^(1/3)); % this root is exactly in the interval [l, u]
        x = real(x) ; % actually, the imaginary part should be zero in theory 
    elseif implement == 3
        % Define the function and its derivative
        f = @(x) x - z - mu ./ (x - l) + mu ./ (u - x);
        df = @(x) 1 + mu ./ ((x - l).^2) + mu ./ ((u - x).^2);
        
        % Newton's method
        x = (l + u) / 2;
        max_iter = 100;
        tol = 1e-12 * length(z);
        for i = 1:max_iter
            x = x - f(x) ./ df(x);
            x = max(x, l + 1e-10);
            x = min(x, u - 1e-10);
            % Check for convergence
            if norm(x - z - mu ./ (x - l) + mu ./ (u - x) ) < tol
                break;
            end
        end
    else
        error('implement should be 1 or 2');
    end

    %% check 
    assert(all(x >= l) && all(x <= u));
    % check by x - z - \mu / ( x - l) + \mu / (u - x) = 0
    % fprintf('check by x - z - mu / ( x - l) + mu / (u - x) = 0, error = %e\n', norm(x - z - mu ./ (x - l) + mu ./ (u - x)));


end


%% on my Macbook Pro , the test result is as follows:
% length(z) = 1e5
% implement 1: takes aound 0.5 seconds
% implement 2: takes around 0.02 seconds
% but the error of implement 2 is larger than implement 1