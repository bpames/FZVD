% % Normalize Coffee and Penicillin data.

x = CoffeeTRAIN(:, 2:287);
[xn, mu, sig] = normalize(x);
mean(xn)
std(xn)

xt = normalize_test(CoffeeTEST(:,2:287),mu,sig);

test = [CoffeeTEST(:,1)+1, xt];
train = [CoffeeTRAIN(:,1)+1, xn];

save('Coffee-normalized.mat', 'train', 'test')

x = CoffeeTRAIN(:, 2:287);
[xn, mu, sig] = normalize(x);
mean(xn)
std(xn)

xt = normalize_test(CoffeeTEST(:,2:287),mu,sig);

test = [CoffeeTEST(:,1)+1, xt];
train = [CoffeeTRAIN(:,1)+1, xn];

save('Coffee-normalized.mat', 'train', 'test')

%% Preprocess penicillium.

clear;
load('penicilliumYES.mat')

% Split into train and test.
trinds = [1:6, 13:18, 25:30];
testinds = [7:12, 19:24, 31:36];

[n,k] = size(Y);
labs = zeros(n,1);
for j=1:n
    for i = 1:k
        if Y(j,i) == 1
            labs(j) = i;
        end
    end
end

Xtr = X(trinds,:);
nzs = any(Xtr,1);
Xtr = Xtr(:, nzs); % Keep only nonzero predictors.
Xte = X(testinds,:);
Xte = Xte(:, nzs); % Keep only nonzero predictors.

labstr = labs(trinds,:);
labste = labs(testinds,:);

%% Normalize X.

[xn, mu, sig] = normalize(Xtr);
mean(xn)
std(xn)

xt = normalize_test(Xte,mu,sig);


%% Copy labels.


    
test = [labste, xt];
train = [labstr, xn];

save('Pen-normalized.mat', 'train', 'test')
