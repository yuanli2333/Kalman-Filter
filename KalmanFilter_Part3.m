%%Project for CCS, part3

%plot the true state position
% figure,
% k=0:1:24;
% plot(10*cos(1/12*k*pi),10*sin(1/12*k*pi),'rs',...
%                        'MarkerEdgeColor','k',...
%                        'MarkerFaceColor','r',...
%                        'MarkerSize',4); 

%录入数据
Y=[7.1165, 9.6022, 8.9144, 9.2717, 6.3400, 4.0484, 0.3411, -0.6784, -5.7726, -5.4925, -9.4523,...
-9.7232, -9.5054, -9.7908, -7.7300, -5.9779, -4.5535, -1.5042, -0.7044, 3.2406,  8.3029, 6.1925, 9.1178, 9.0904, 9.0662];
Z=[0.000, 3.1398,  6.3739, 9.5877, 10.1450, 10.1919,9.0683, 10.2254, 7.5799, 7.7231,  5.4721,...
3.3990,  0.9172, -1.3551,  -5.2708, -9.7011, -9.4256, -9.3053, -9.3815, -9.8822, -8.1876,-8.7501,-4.5653, -1.9179, -1.000];
%得到测量值
y=[Y;Z];
I_25=ones(25,1);

%最小二乘法求半径和圆心
Q=[2*Y',2*Z',I_25];
T=[Y.^2+Z.^2]';
X_shu=inv(Q'*Q)*Q'*T;
R = sqrt(X_shu(3));
%最小二乘法求速度theta
Q_t = linspace(0,24,25)';
T_t = [mod(acos(Y/R)',2*pi) mod(asin(Z/R)',2*pi)];
theta = (Q_t'*Q_t)^(-1)*Q_t'*T_t; 

%%Kalman Filter
N=25;
Nrum=1;
sw=2;
sv=1;
%A=[cos(theta(2)) -sin(theta(2)); sin(theta(2)) cos(theta(2))];
A=[1 0;0 1];
C=[1 0;0 1];
R1=sw^2*[1 0;0 1];
%R1=[0.3290 0.2773]'*[0.3290 0.2773];
%R1=[0 0;0 0];
R2=sv^2*[1 0;0 1];
Pm(:,:,1)=1e7*eye(2); %P(0|-1),Pm is the P(k/k-1), P is the P(k/k)

%Caculation of Kalman Gains
for k=1:N
    Kf(:,:,k)=Pm(:,:,k)*C'*(C*Pm(:,:,k)*C'+R2)^(-1); %Kf(:,k)=Kf(k)
    K(:,:,k)=(A*Pm(:,:,k)*C')*(C*Pm(:,:,k)*C'+R2)^(-1); %K(:,k)=K(k)
    P(:,:,k)=Pm(:,:,k)-(Pm(:,:,k)*C')*(C*Pm(:,:,k)*C'+R2)^(-1)*C*Pm(:,:,k); %P(:,:,k)=P(N|N)
    Pm(:,:,k+1)=A*Pm(:,:,k)*A'+R1-K(:,:,k)*(C*Pm(:,:,k)*C'+R2)*K(:,:,k)'; %Pm(:,:,k+1)=P(N+1|N)
end


%Innitialization
w1=random('normal',0,sw,Nrum,N); %process noise
w2=random('normal',0,sw,Nrum,N); %process noise
v1=random('normal',0,sv,Nrum,N); %measurement noise
v2=random('normal',0,sv,Nrum,N); %measurement noise

xhm(:,1)=[10 0]'; %xhm(:,1)=x(0|-1)
x(:,1)=[10 0]'; %x(:,1)=x(0)   


for k=1:N
    %the model
    x(:,k+1)=A*x(:,k)+[w1(k);w2(k)];
    %y(:,k)是测量值，直接由上面的 y=[Y;Z]可得，无需使用模型计算;
%     y(:,k)=C*x(:,k)+[v1(k);v2(k)];
    %Kalman Filter
    %if k==1
        xh(:,k)=xhm(:,k)+Kf(:,:,k)*(y(:,k)-C*xhm(:,k)); %xh(:,k)=x^hat(k|k)
    %else
    %    xh(:,k)=(1-Kf(:,:,k))*A*xh(:,k-1)+Kf(:,:,k)*y(:,k);
   % end
    
    xhm(:,k+1)=A*xhm(:,k)+K(:,:,k)*(y(:,k)-C*xhm(:,k)); %xhm(:,k+1)=x^hat(k+1|k)
end
  

%画图
figure,
k=0:1:N-1;
h1=plot(10*cos(1/12*k*pi),10*sin(1/12*k*pi),'rs',...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',4); 
hold on,
h2=plot(Y,Z,'.b','LineWidth',1);
h3=plot(xh(1,:),xh(2,:),'om:');
legend([h1(1),h2(1),h3(1)], 'True state', 'Measurements','KalmanFilter');
axis equal;
xlabel('y-axis');
ylabel('z-axis');
title('Part3. By YuanLi');

% caculate Biases and variances
%先计算出true state position
k=0:(N-1);
Tr_y=10*cos(2*pi*k/(N-1)); % true state positions for y 
Tr_z=10*sin(2*pi*k/(N-1)); % true state positions for x

Bias_y=1/(N-1)*sum(Tr_y - xh(1,:));
Bias_z=1/(N-1)*sum(Tr_z - xh(2,:));

Var_y=(1/(N-1))*sum((Tr_y - xh(1,:)).^2);
Var_z=(1/(N-1))*sum((Tr_z - xh(2,:)).^2);

Bias_y
Bias_z
Var_y
Var_z

