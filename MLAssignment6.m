clear all
close all
clc

%qq = csvread('house_prices_data_training_data.csv',1,2);

data = csvread('house_prices_data_training_data.csv',1,2);

Y = data(:,1);
data = data(:,2:end);

Corr_x = corr(data);
x_cov = cov(data) ;
[U S V] =  svd(x_cov);
eigenValues = diag(S);
m = length(data(1,:));



for k=1:1:m
    numerator = sum(eigenValues(1:k));
    denimerator = sum(eigenValues);
    alpha = 1 - (numerator / denimerator)
    if (alpha <= 0.001)
        ss = k;
        break
    end
end

Reduced_Data =  U(:,1:k)' * data';
approximateData = U(:,1:k) * Reduced_Data;

n = length(approximateData(1,:));
error = (1/m) * (data - approximateData');

%Y=(Y-mean(Y))/std(Y);
[ iterations , Js ,theta] = helperMLAss1(Reduced_Data' , alpha , Y);











% for w=2:n
%     if max(abs(approximateData(:,w)))~=0
%     approximateData(:,w)=(approximateData(:,w)-mean((approximateData(:,w))))./std(approximateData(:,w));
%     end
% end
%
% Y=data(:,3)/mean(data(:,3));
%
% Theta=zeros(n,1);
% k=1;
%
% E(k)=(1/(2*m))*sum((approximateData*Theta-Y).^2);
%
% R=1;
% while R==1
% alpha=alpha*1;
% Theta=Theta-(alpha/m)*approximateData'*(approximateData*Theta-Y);
% k=k+1
% E(k)=(1/(2*m))*sum((approximateData*Theta-Y).^2);
% if E(k-1)-E(k)<0
%     break
% end
% q=(E(k-1)-E(k))./E(k-1);
% if q <.000001;
%     R=0;
% end
% end





%part 2
Reduced_Data = Reduced_Data';
for randomk=2:3
    
    
    %randomk = 4
    
    for i=2:3
        for j=1:1:randomk
            centroids(j,i) = Reduced_Data(ceil(21607*rand(1,1)) , ceil(3*rand(1,1)));
        end
    end
    
    
    normalisedData = normalizedata(Reduced_Data,3);
    for iterations=1:25
        
        for i=1:1:21607
            for j=1:1:randomk
                error2(j) = sum(sum((normalisedData(i,:) - centroids(j,:)).^2));
            end
            minError2(randomk) = min(error2);
            index(i) = find(error2 == minError2(randomk),1);
        end
        
        
        for i=1:1:randomk
            indexes =  find(index(i) == i);
            if(size(indexes,1) > 0)
                centroids(i,:) = mean(normalisedData(indexes,i))
            end
        end
        
    end
    
end

figure(1)
plot(1:randomk,minError2)








%%
for randomk=2:10
    
    
    %randomk = 4
    
    for i=1:1:18
        for j=1:1:randomk
            centroids(j,i) = data(ceil(21607*rand(1,1)) , ceil(18*rand(1,1)));
        end
    end
    
    
    normalisedData = normalizedata(data,18);
    for iterations=1:50
        
        for i=1:1:21607
            for j=1:1:randomk
                error2(j) = sum(sum((normalisedData(i,:) - centroids(j,:)).^2));
            end
            minError2(randomk) = min(error2);
            index(i) = find(error2 == minError2(randomk),1);
            dataError(i) = minError2(randomk);
        end
        
       resultError(randomk) = sum(dataError);
        for i=1:1:randomk
            indexes =  find(index(i) == i);
            if(size(indexes,1) > 0)
                centroids(i,:) = mean(normalisedData(indexes,i))
            end
        end
        
    end
    
end
figure(2)
plot(2:randomk,resultError(2:randomk))



%%
%part 3

meann = mean(data);
standardDeviation = std(data);
counterPart3 = 0;

for q=1:21607
    
    for i=1:1:18
        pdf(i) = normcdf(data(q,i),mean(i),standardDeviation(i));
    end
    
    pdfProduct = prod(pdf(1:5));
    pdfProduct = pdfProduct * prod(pdf(7:end));
    
    status = 'Non Anomaly';
    
    if(pdfProduct > 0.999) || (pdfProduct < 10^(-50))
        status = 'Anomaly';
        counterPart3 = counterPart3 + 1;
    end
    
end
