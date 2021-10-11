clear all;
%WSN parameters:
N = 35;   %the number of slave nodes;
K = 10;   %the number of the tasks in each DAG; 

Bati = randi([1,N],N,1);  % the battery energy of slave nodes;
Batm = randi([N,2*N],1,1); % the battery energy of master nodes;

Ni = zeros(N);             % Neighbor nodes set;
Ni(1,2) = 1; Ni(1,3) = 1;Ni(1,5) = 1;
Ni(2,1) = 1; Ni(2,4) = 1;
Ni(3,1) = 1; Ni(3,5) = 1;
Ni(4,2) = 1;
Ni(5,1) = 1; Ni(5,3) = 1;

% The number of tasks in one app;
Num_Task = [5,10,15,20];
Taskselect = [5,10,16,21];
% Matrix represent the DAG flow;
ori_DAG = [-1,1,zeros(1,18);
    0,-1,1,zeros(1,17);
    0,0,-1,1,zeros(1,16);
    0,0,-1,0,1,zeros(1,15);
    0,0,0,-1,1,zeros(1,15);
    0,0,0,0,-1,1,zeros(1,14);
    0,0,0,0,0,-1,1,zeros(1,13);
    0,0,0,0,0,0,-1,1,zeros(1,12);
    zeros(1,7),-1,1,zeros(1,11);
    zeros(1,8),-1,1,zeros(1,10);
    zeros(1,9),-1,1,zeros(1,9);
    zeros(1,9),-1,0,1,zeros(1,8)
    zeros(1,10),-1,1,zeros(1,8);
    zeros(1,11),-1,1,zeros(1,7);
    zeros(1,12),-1,1,zeros(1,6);
    zeros(1,13),-1,1,zeros(1,5);
    zeros(1,14),-1,1,zeros(1,4);
    zeros(1,15),-1,1,zeros(1,3);
    zeros(1,16),-1,1,zeros(1,2);
    zeros(1,17),-1,1,0;
    zeros(1,18),-1,1];

DAG = ori_DAG(1:Taskselect(2),1:Num_Task(2));
    
save wsndata N K Ni Bati Batm DAG;