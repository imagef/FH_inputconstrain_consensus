function [K_bar,K,A,B,N]=fh_consensus_constraim(eps)

%系统方程
% A = [0 1 0 0;
%     0 0 1 0;
%     0 0 0 1;
%     -1 0 -2 0];
A = [0 1;-1 1];
% A=[0.995,0.09983;-0.09983,0.995];
% A=[0,1;0.1,1.05];
% A=[0,1;-1,0.5];

% E = [0.0621; -0.6784;0.1655;0.4247];
% B = [-1;1;1;3];
% A = [1 1;0 1];
% A = [0.3 0.16; 0.92 0.82];
% B = [-1;0];

% A = [1 1; 0 1];
B = [-1;0];
% B = [0;1];




R = 1*eye(1);

N =101;
Step = 100;

% eps = 2;
Q = eps*diag([1 1], 0);

for i = 1:Step 
    u = rand(1,N-1);
    x(1:2,i,1) = randn(2,1);
    for k = 1:N-1
       x(1:2,i,k+1) = A * x(1:2,i,k) + B * u(:,k) ; % 状态量
       x(3,i,k) = u(:,k);
    end
end

for k = 1: N-1
    for i = 1: Step 
        fai(k,i,:) = [x(1,i,k)^2  2*x(1,i,k)*x(2,i,k) 2*x(1,i,k)*x(3,i,k) x(2,i,k)^2 2*x(2,i,k)*x(3,i,k) x(3,i,k)^2];
    end
end

for k = N-1:-1:1 

    if k == N-1   
        for i = 1:Step
            gamma(k,i) = x(1:2,i,k+1)' * Q* x(1:2,i,k+1) + x(1:2,i,k)' * Q * x(1:2,i,k) + x(3,i,k)' * R * x(3,i,k);  
        end
    else
        for i = 1:Step
            % gamma(k,i) = [x(1:2,i,k+1); 4/30*K_bar(:,k+1)'*x(1:2,i,k+1)]' * H(:,:,k+1) * [x(1:2,i,k+1); 4/30*K_bar(:,k+1)'*x(1:2,i,k+1)] + x(1:2,i,k)' * Q * x(1:2,i,k) + x(2,i,k)' * R * x(2,i,k);
            gamma(k,i) = x(1:2,i,k+1)' * P_bar(:,:,k+1) * x(1:2,i,k+1) + x(1:2,i,k)' * Q * x(1:2,i,k) + x(3,i,k)' * R * x(3,i,k);
        end
    end
    rank(squeeze(fai(k,:,:))'*squeeze(fai(k,:,:)))
    lamda(:,k) = inv(squeeze(fai(k,:,:))'*squeeze(fai(k,:,:))) * squeeze(fai(k,:,:))' * gamma(k,:)';
    H_xx(:,:,k) = [lamda(1,k) lamda(2,k); lamda(2,k) lamda(4,k)];
    H_xu(:,:,k) = [lamda(3,k);lamda(5,k)];
    H_ux(:,:,k) = H_xu(:,:,k)';

    H_uu(k) = lamda(6,k);

    H(:,:,k) = [H_xx(:,:,k) H_xu(:,:,k);  H_ux(:,:,k) H_uu(k)];

    K_bar(:,k) = -inv(H_uu(k))*H_ux(:,:,k);
    P_bar(:,:,k) = H_xx(:,:,k) + 4/30*K_bar(:,k)*H_ux(:,:,k)+H_xu(:,:,k)*4/30*K_bar(:,k)'+4/30*K_bar(:,k)*H_uu(k)*4/30*K_bar(:,k)';
    NormFP_bar(k) = norm(P_bar(:,:,k),'fro');
    
end


for kk = N-1:-1:1
    if kk == N-1
        P(:,:,kk) = Q+A'*Q*A-A'*Q*B*inv(R+B'*Q*B)*B'*Q*A;
    else
        P(:,:,kk) = Q+A'*P(:,:,kk+1)*A - (2*4/(6*5)-(4/(6*5))^2)*A'*P(:,:,kk+1)*B*inv(R+B'*P(:,:,kk+1)*B)*B'*P(:,:,kk+1)*A;
    end
    NormFP(kk) = norm(P(:,:,kk),'fro');
end

for kk = 1:N-1
    if kk == N-1
        K(:,kk) = inv(R + B'*Q*B)*B'*Q*A; %应该是K(N-1)
    else
        K(:,kk) = inv(R + B'*P(:,:,kk+1)*B)*B'*P(:,:,kk+1)*A; % 这里应该是P(k+1)，但是
    end
end

