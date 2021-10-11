% clear all;
%WSN parameters:
N = 3;   %the number of slave nodes;
K = 5;   %the number of the tasks in each DAG;
var = 1; %the variation among the tasks;
Bati = [2,5,1];                    % randi([1,3],N,1);
Batm = 6;                                                              

%Processing Energy Cost Parameters:
P_i = 36.9e-3; P_m = 36.9e-3;          % Processing Power; unit: W
f_i = 32e6; f_m = 32e6;                % Processing speed; unit: Hz

%Communication Energy Cost Parameters;
P_T0 = 59.8e-3;                        % RF circuit energy cost for receiving and transmitting; unit: W
t_tx = 4e-6; t_rx = 4e-6;              % time for transmitting or receiving 1 bitt data packet; unit: s

epi = P_i*(1/f_i);                % parameters of processing energy without workload of slave nodes;
epm = P_m*(1/f_m);                % parameters of processing energy without workload of master nodes;

eci = P_T0*(t_tx);                     % parameters of communication energy without workload of slave nodes;
ecm = P_T0*(t_rx);                     % parameters of communication energy without workload of master nodes;

%computation workloads of each DAG graph;
Wk = randi([100,1000],N,K);
% Wk = [415,655,848,926,779;277,526,627,357,442;326,416,595,782,611];
%the amount of the communicated data on the edges;
le = randi([100,500],N,K);
% le = [130,412,328,235,224;121,474,288,165,311;312,152,104,418,166];

% Wk = [516,827,632,978,523;
%       835,593,905,516,631;
%       768,614,598,943,601];
% le = [500,300,200,400,500;
%       500,300,200,400,500;
%       500,300,200,400,500;];
% Wk = [834,922,350,969,962;916,669,592,242,537;214,187,962,974,821];
% le = [156,417,114,372,257;269,484,440,403,362;467,362,474,397,168];
% Wk = [736,141,726,131,789;128,187,385,495,816;349,841,956,443,268];
% le = [296,384,372,147,236;278,402,362,299,334;359,210,165,484,189];
Wk = [196,801,902,278,550;689,744,401,127,532;545,914,729,770,915];
le = [344,423,196,296,385;347,331,455,167,300;444,173,111,492,288];
%incidence matrix of the graph:
DAG = [-1,1,0,0,0;
        0,-1,0,1,0;
        0,-1,1,0,0;
        0,0,-1,1,0;
        0,0,0,-1,1];
[row,column] = size(DAG);
%Ln
Ln = le*(-DAG);

%% DAG of sum of xii and two of zji;
dag_sum = kron(ones(1,N),DAG);
Dag_sum = kron(eye(N),dag_sum);
DAG_sum = [zeros(N*row,1) Dag_sum zeros(N*row,K*N*(N-1))];

% B*Xii;
dag_x = zeros(N,N*N);
for i = 1:N
    dag_x(i,N*(i-1)+i) = 1;
end
Dag_x = kron(dag_x,DAG); 

% constriant: DAG of xii;
DAG_x = [zeros(N*row,1) Dag_x zeros(N*row,K*N*(N-1))];

% constriant: DAG of sum of xii and zji;
Dag = matrixstructure(N);
xzij = kron(Dag,DAG);
xzij = [zeros(N*(N-1)*row,1),xzij,zeros(N*(N-1)*row,K*N*(N-1))];
Constriant1 = [DAG_sum;DAG_x; xzij];

%constriant:sum of x and z are smaller than 1;
sum_xz = [zeros(N*(N-1)*K,1) kron(Dag,eye(K)) zeros(N*(N-1)*K,K*N*(N-1))];
B0 = ones(N*(N-1)*K,1);

Constriant1_1 = kron(eye(N-1),DAG);
Constriant1_2 = kron(eye(N),Constriant1_1);
Constriant1_3 = [zeros(N*(N-1)*row,N*N*K+1),Constriant1_2];

Constriant1 = [Constriant1;Constriant1_3];
B1 = zeros(2*N*row+2*(N*(N-1)*row),1);

%B*Zji and Zjii;

Zji =[zeros(N,1) eye(N)];
ZJI = [kron(eye(N-1),Zji) zeros(N*(N-1),1)];
Zjii = eye(N*(N-1));
sumzz = [ZJI Zjii];
SUMzz = [zeros(N*(N-1)*row,1) kron(sumzz,DAG)];
bsum = zeros(N*(N-1)*row,1);

% sum(Zji)-Zjii(1) >=0;

sumzji = kron(ZJI, -1*ones(1,K));
zjii_1 = kron(Zjii,[1,zeros(1,K-1)]);
sum_1 = [zeros(N*(N-1),1) sumzji zjii_1];
bsum_1 = zeros(N*(N-1),1);

% The sum of Zji and Zjii is less than one;
sumz = [zeros(N*(N-1)*K,1) kron(sumzz,eye(K))];
bsumz = ones(N*(N-1)*K,1);

%% T>=Ei/Bati;

%parameter of Xii:
Xii = epi*Wk*10^3 + eci*Ln;
Constriant2 = [];
for i = 1:N
    repeat = kron(eye(N),Xii(i,:));
    Constriant2 = [Constriant2 repeat];
end

%parameter of Zijj;
Zijj = 2*eci*Ln;
Constriant3 = [];
for i = 1:N
    repeat1 = kron(eye(N),Zijj(i,:));
    repeat1(:,(i*K-K+1):i*K) = [];
    Constriant3 = [Constriant3 repeat1];
end

Constriant4 = [Constriant2 Constriant3];
for i = 1:N
    Constriant4(i,:)=Constriant4(i,:)./Bati(i);
end
 
Constriant4 = [-1*ones(N,1) Constriant4];
B2 = zeros(N,1);

%% ===================================================================
Zi=[];
for i = 1:N
    Z = zeros(N*K,N*K);
    Z((i*K-K+1):i*K,(i*K-K+1):i*K) = -1*eye(K);
    Zi = [Zi Z];
end

Zj = kron(ones(1,N-1),eye(K));
Zj_new = kron(eye(N),Zj);

Constriant5 = [zeros(N*K,1) Zi Zj_new];
B3 = zeros(N*K,1);
%% =============================================================
% Em;
Em = ecm*Ln - epm*Wk*10^3;
Emt = kron(Em,ones(N,1));
Emt = reshape(Emt',1,[]);
Constriant6 = [Emt,zeros(1,K*N*(N-1))]./Batm;
Constriant6 = [-1 Constriant6];

Epm = epm*Wk*10^3;
B4 = -sum(sum(Epm))/Batm;


%%======================================================== 
% Equation
Node1 = [1,zeros(1,K-1)];

Node_first =[];
for i = 1:N
    tool = zeros(N);
    tool(i,i) = 1;
    Node1_new = kron(tool,Node1);
    Node_first = [Node_first Node1_new];
end
equation1 = [zeros(N,1) Node_first zeros(N,K*N*(N-1))];
B_equal1 = ones(N,1);

Node0 = [zeros(1,K-1),1];
Node0_new = kron(eye(N),Node0);
Node0_new = kron(eye(N),Node0_new);

% for i = 1:N
%     Node0_new(1+N*(i-1),:) =[];
% end

equation2 = [zeros(N*N,1) Node0_new zeros(N*N,K*N*(N-1))];
B_equal2 = zeros(N*N,1);



%%================================
lb = zeros(K*N*(2*(N-1)+1)+1,1);
ub = [inf;ones(K*N*(2*(N-1)+1),1)];

%%====================================
% objective function

f = [1;zeros(K*N*(2*(N-1)+1),1)];
A = [sum_1;sumz;SUMzz;sum_xz;Constriant1;
    Constriant4;
    Constriant5;
    Constriant6];
b = [bsum_1;bsumz;bsum;B0;B1;B2;B3;B4];
intcon = (2:K*N*(2*(N-1)+1)+1);
Aeq = [equation1;equation2];
beq = [B_equal1;B_equal2];

[x1,fval1] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
x1(1)=[];
z = x1(N*N*K+1:K*N*(2*(N-1)+1));
X = x1(1:N*N*K);
X = reshape(X,N*K,N);
X = X';

Constriant4(:,1) = [];
Constriant6(1) = [];

Ei = Constriant4 * x1;
EM = Constriant6 * x1 - B4;

%% Task allocation without algorithm

q = randi([1,K-1]);
xnor = [ones(q,1);zeros(K-q,1)];
Xnor = kron(ones(1,N),xnor);
E_nor =[];
em = [];
for i = 1:N
    E_nor = [E_nor,Xii(i,:)*Xnor(:,i)/Bati(i)];
    em = [em, Em(i,:)*Xnor(:,i)/Batm];
end
Em = sum(em)-B4;
E_nor = [E_nor, Em];
Emin = max(E_nor);

increase = Emin/fval;
% %================================================================================================
% %% verify 
% A_ij = [1,0;1,0;1,0];
% [x_final_1] = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;1,0;0,0];
% x_final_2 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;1,0;0,1];
% x_final_3 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,0;1,0];
% x_final_4 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,0;0,0];
% x_final_5 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,0;0,1];
% x_final_6 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,1;1,0];
% x_final_7 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,1;0,0];
% x_final_8 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [1,0;0,1;0,1];
% x_final_9 = multihopverify(Wk,le,DAG,A_ij);
% 
% A_ij = [0,0;1,0;1,0];
% x_final_10 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;1,0;0,0];
% x_final_11 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;1,0;0,1];
% x_final_12 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,0;1,0];
% x_final_13 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,0;0,0];
% x_final_14 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,0;0,1];
% x_final_15 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,1;1,0];
% x_final_16 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,1;0,0];
% x_final_17 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,0;0,1;0,1];
% x_final_18 = multihopverify(Wk,le,DAG,A_ij);
% 
% A_ij = [0,1;1,0;1,0];
% x_final_19 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;1,0;0,0];
% x_final_20 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;1,0;0,1];
% x_final_21 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,0;1,0];
% x_final_22 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,0;0,0];
% x_final_23 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,0;0,1];
% x_final_24 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,1;1,0];
% x_final_25 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,1;0,0];
% x_final_26 = multihopverify(Wk,le,DAG,A_ij);
% A_ij = [0,1;0,1;0,1];
% x_final_27 = multihopverify(Wk,le,DAG,A_ij);
% 
% x_final = [x_final_1 x_final_2 x_final_3 x_final_4 x_final_5 x_final_6 x_final_7 x_final_8 x_final_9 x_final_10 x_final_11 x_final_12 x_final_13 x_final_14 x_final_15 x_final_16 x_final_17 x_final_18 x_final_19 x_final_20 x_final_21 x_final_22 x_final_23 x_final_24 x_final_25 x_final_26 x_final_27];
% [min_value,location] = min(x_final(1,:));

