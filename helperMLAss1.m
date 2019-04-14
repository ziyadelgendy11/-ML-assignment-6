function [ iterations , Js ,theta] = helperMLAss1( X , alpha , results)

n = size(X(1,:));
theta = zeros(n(2),1);
%Normalize all the features
for w=2:n(2)
    if max(abs(X(:,w)))~=0
    X(:,w)=(X(:,w)-mean((X(:,w))))./std(X(:,w));
    end
end


%Normalize the result (prices)
results(:)=(results(:)-mean((results(:))))./std(results(:));


% Compute the Cost Function
%X = [ones(length, 1), feature1 , feature2, feature3];
J = ComputeCost(X, results, theta);


% Run Gradient Descent
[theta, Js] = GradientDescent(X, results, theta, alpha);


% plotting the cost function
iterations=size(Js);


end

