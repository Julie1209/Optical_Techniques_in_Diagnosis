x_fivetimes = zeros(5, 10000);
N = zeros(5, 20);

figure;

for i = 1:5
    x = rand(1, 10000);
    nbins = 20;
    edges = 0:0.05:1;
    Y = discretize(x, edges);
    N(i, :) = histcounts(Y);
    
    subplot(2, 3, i);
    histogram(x, nbins);
    title(['Run ', num2str(i)]);
    xlabel('Interval');
    ylabel('Number of random numbers');
end
x_mean = mean(N);
disp(x_mean);
x_std = std(N);
disp(x_std);

