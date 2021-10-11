%% The number of devices;
N = [5,10,15,20,25,30,35,40];
N_Y = load('result_N');
bar(N,Y);
title('Effect of the number of end device in application');
xlabel('Number of End Device,N');
ylabel('Network Lifetime Increase')

%% The number of tasks in one application;
K = [5,10,15,20];
K_Y = load('result_K');
bar(K_Y,Y);
title('Effect of the number of tasks in application');
xlabel('Number of Tasks,K');
ylabel('Network Lifetime Increase')

%% The variation level in one application;
var = {'low','middle','high'};
var_Y = load('result_var');
bar(var_Y);
set(gca, 'XTickLabel', var);
title('Effect of the variation of tasks in application');
xlabel('Variation levels among the tasks');
ylabel('Network Lifetime Increase')