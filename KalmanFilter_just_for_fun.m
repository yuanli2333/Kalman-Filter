%% just for fun

clear
N=100; %number of samples, if N很大，到一定阶段，Kf,K,P都变成constant
Nrun=1; %number of runs
T=1; %sampling interval,采样间隔
sw=0; % sigma_w
sv=1; % sigma_v
%%%%%%%%%
% Process
%%%%%%%%%
A=[1 T; 0 1]; %process matrix
C=[1 0]; %process matrix
R1=sw^2*[T^4/4 T^3/2; T^3/2 T^2]; %covariance of process noise
R2=sv^2; %covariance of measurement noise
Pm(:,:,1)=1e5*eye(2); %P(0|-1),Pm is the P(k/k-1), P is the P(k/k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Kalman Gains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
    Kf(:,k)=Pm(:,:,k)*C'*(C*Pm(:,:,k)*C'+R2)^(-1); %Kf(:,k)=Kf(k)
    K(:,k)=(A*Pm(:,:,k)*C')*(C*Pm(:,:,k)*C'+R2)^(-1); %K(:,k)=K(k)
    P(:,:,k)=Pm(:,:,k)-(Pm(:,:,k)*C')*(C*Pm(:,:,k)*C'+R2)^(-1)*C*Pm(:,:,k); %P(:,:,k)=P(N|N)
    Pm(:,:,k+1)=A*Pm(:,:,k)*A'+R1-K(:,k)*(C*Pm(:,:,k)*C'+R2)*K(:,k)'; %Pm(:,:,k+1)=P(N+1|N)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializaton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=random('normal',0,sw,Nrun,N); %process noise
v=random('normal',0,sv,Nrun,N); %measurement noise
xhm(:,1)=[0 0]'; %xhm(:,1)=x(0|-1)
x(:,1)=[0 1]'; %x(:,1)=x(0)

for run=1:Nrun
    for k=1:N
        %%%%%%%%%
        % Process
        %%%%%%%%%
        x(:,k+1)=A*x(:,k)+[T^2/2 T]'*w(run,k); %x(:,k+1)=x(k+1)
        y(k)=C*x(:,k)+v(run,k); %y(k)=y(k)
        %xh(:,0)=xhm(:,1)+Kf(0)*(y(0)-C*xhm(:,1));
        %%%%%%%%%%%%%%%
        % Kalman Filter
        %%%%%%%%%%%%%%%
        xh(:,k)=xhm(:,k)+Kf(:,k)*(y(k)-C*xhm(:,k)); %xh(:,k)=x^hat(k|k)
        %xh(:,k)=(1-Kf(:,k))*a*xh(:,k-1)+Kf(:,k)*y(k);
        xhm(:,k+1)=A*xhm(:,k)+K(:,k)*(y(k)-C*xhm(:,k)); %xhm(:,k+1)=x^hat(k+1|k)
    end
end

% ax=1:N-1;
% plot(ax,x(1,1,N));


%%for moive
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated1.gif';

for i=2:40
    plot(x(1,i),10,'-ro','LineWidth',1,...
                           'MarkerEdgeColor','k',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',12);
    xlabel('k');
    axis equal;
    axis([0 40 0 20]);
    title('Just for fun: Life just like KalmanFilter. By YuanLi');
    hold all
    plot(x(1,i-1),10,'-g*','LineWidth',1,...
                           'MarkerEdgeColor','k',...
                           'MarkerFaceColor','g',...
                           'MarkerSize',20);
    hold all
    plot(y(i),10,':bs','LineWidth',1,...
                           'MarkerEdgeColor','k',...
                           'MarkerSize',20);  
    %gtext('solid-line: x(k)');
    hold off
    %M(i-1)=getframe;
   frame = getframe(h); 
   im = frame2im(frame); 
   [imind,cm] = rgb2ind(im,256); 
   % Write to the GIF File  
end
imwrite(imind,cm,filename,'gif','Loopcount',inf)
%movie(M,1,3);
%movie2avi(M,'KalmanFilter.avi');


