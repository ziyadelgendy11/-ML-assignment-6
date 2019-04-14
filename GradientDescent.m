function [theta, Js] = GradientDescent(X, y, theta, alpha)

% Prepare Variables
    m = length(y);
    %Js = zeros(iterations, 1);
    checker = true;
    counter = 1;
    
    while checker == true 
        
        h = X * theta;
        
       %theta(1) = theta(1) - (alpha * (1 / m) * sum(h - y));
       % for j = 2 : length(theta)
           % theta(j) = theta(j) - (alpha * (1 / m) * sum((h - y) .* X(:, j)));
       % end
        
       theta=theta-(alpha/m)*X'*(X*theta-y);
       Js(counter) = ComputeCost(X, y, theta);
        
        
        if(counter>1)
            if( (Js(counter-1)-Js(counter)) / Js(counter-1) < 0.000001)
                checker = false;
            end
            
           % if Js(counter-1)-Js(counter)<0
           %     checker = false;
           % end 
            
        end
        counter = counter + 1;
    end

end

