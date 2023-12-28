clear all;
clc;
format long;

% Assumes constant velocity linear motion
% set real position that linear motion along (y=x) line at time interval
realPOS = zeros(2,11);
for j = 1:1:11
    realPOS(1,j) = j-1; realPOS(2,j) = j-1;
end

% define struct
ToA = struct;
ToA.position = struct;
ToA.position.real = realPOS;

% set anchor position
ToA.position.anchor = [0 10; 0, 0; 10, 0; 10, 10];

% set estimation position by ToA situation simulation
ToA.position.hat = zeros(2,11);

% set estimation position by Kalman filter with ToA situation simulation
ToA.aPos_state_s = zeros(10,11);
ToA.aPri_state_s = zeros(10,11);
ToA.my_aPri_state_s  = zeros(10,11);
ToA.my_aPos_state_s  = zeros(10,11);

% set errors
ToA.error_s = zeros(5,11);
ToA.error_v_k_s = zeros(5,11);
ToA.error_kf_s = zeros(5,11);
ToA.my_error_kf_s = zeros(4,11);

% noise with a normal distribution with a variance of sigma^2
sigma_index = [0.01 0.1 1 10 100];

% number of simulation iterations
iter = 1e4;

% Optimal Q, w values obtained by repeating simulations
Q_s_and_w_s;

% time interval
deltaTime = 0.1;

% in ToA ... H * [state vector]' = Z (Z is equal to a distance term) 
% set coefficient matrix H in ToA
H = ...
   [-2*(ToA.position.anchor(1,:)-ToA.position.anchor(2,:));
    -2*(ToA.position.anchor(1,:)-ToA.position.anchor(3,:));
    -2*(ToA.position.anchor(1,:)-ToA.position.anchor(4,:));
    -2*(ToA.position.anchor(2,:)-ToA.position.anchor(3,:));
    -2*(ToA.position.anchor(2,:)-ToA.position.anchor(4,:));
    -2*(ToA.position.anchor(3,:)-ToA.position.anchor(4,:))
    ];

% set prediction coefficient matirx of kalman filter as [1 0; 0 1]
A = eye(1,1);

for s = 1:1:5
    % set parameters 
    sigma = sigma_index(s);

    %Set up state for iteration
    aPos_state = zeros(2,11);
    aPri_state = zeros(2,11);
    my_aPos_state = zeros(2,11);
    my_aPri_state = zeros(2,11);

    %Set up errors for iteration
    error = zeros(1,11);
    error_v_k = zeros(1,11);
    error_kf = zeros(1,11);
    my_error_kf = zeros(1,11);
   
    %% iteration part of a sigma(i)
    % set parameters

    mean_w = mean(w_s_ram(2*s-1:2*s,4:11),2);
    Q_k = cell2mat(Q_s_ram(s,1));
    
    %iteration
    for i = 1:1:iter

        %Set up state for iteration
        hatPOS = zeros(2,11);
        hat_v_k_POS = zeros(2,11);

        aPos_state_k=zeros(2,11);
        aPri_state_k=zeros(2,11);     
        my_aPos_state_k=zeros(2,11);
        my_aPri_state_k=zeros(2,11);

        %Set up errors for iteration        
        error_k = zeros(1,11);
        error_v_k_k = zeros(1,11);
        error_kf_k = zeros(1,11);
        my_error_kf_k = zeros(1,11);

        % ToA tracking modeling and simulation
        for j = 1:1:11
            noise = zeros(1,4); distance = zeros(1,4); distance_hat = zeros(1,4);
            
            % set of a WGN in ToA, a hat distance for a sigma(i)
            for k = 1:1:4
                noise(k) = (sigma^2) * randn;
                distance(k) = norm(ToA.position.real(:,j)-ToA.position.anchor(k,:)');
                distance_hat(k) = noise(k) + distance(k);
            end
            
            % ToA tracking
            Z = ...
           [(distance_hat(1))^2-(distance_hat(2))^2-ToA.position.anchor(1,1)^2+ToA.position.anchor(2,1)^2-ToA.position.anchor(1,2)^2+ToA.position.anchor(2,2)^2;
            (distance_hat(1))^2-(distance_hat(3))^2-ToA.position.anchor(1,1)^2+ToA.position.anchor(3,1)^2-ToA.position.anchor(1,2)^2+ToA.position.anchor(3,2)^2;
            (distance_hat(1))^2-(distance_hat(4))^2-ToA.position.anchor(1,1)^2+ToA.position.anchor(4,1)^2-ToA.position.anchor(1,2)^2+ToA.position.anchor(4,2)^2;
            (distance_hat(2))^2-(distance_hat(3))^2-ToA.position.anchor(2,1)^2+ToA.position.anchor(3,1)^2-ToA.position.anchor(2,2)^2+ToA.position.anchor(3,2)^2;
            (distance_hat(2))^2-(distance_hat(4))^2-ToA.position.anchor(2,1)^2+ToA.position.anchor(4,1)^2-ToA.position.anchor(2,2)^2+ToA.position.anchor(4,2)^2;
            (distance_hat(3))^2-(distance_hat(4))^2-ToA.position.anchor(3,1)^2+ToA.position.anchor(4,1)^2-ToA.position.anchor(3,2)^2+ToA.position.anchor(4,2)^2];    
            
            hatPOS(:,j) = (inv(H' * H) * H') * Z;

            % part of a ToA Tracking with KF 
            v_k =...
           [-2*(distance_hat(1)*noise(1))+2*(distance_hat(2)*noise(2));
            -2*(distance_hat(1)*noise(1))+2*(distance_hat(3)*noise(3));
            -2*(distance_hat(1)*noise(1))+2*(distance_hat(4)*noise(4));
            -2*(distance_hat(2)*noise(2))+2*(distance_hat(3)*noise(3));
            -2*(distance_hat(2)*noise(2))+2*(distance_hat(4)*noise(4));
            -2*(distance_hat(3)*noise(3))+2*(distance_hat(4)*noise(4))];       
            
            hat_v_k_POS(:,j) = (inv(H' * H) * H') * (Z+v_k);
            Z_KF = H * hat_v_k_POS(:,j);
            
            %% with KalmanFilter 
            % prediction part !!
            if j == 3
                aPos_state_k(:,j) = hatPOS(:,j) ;   
                aPos_P = (realPOS(:,j)-hatPOS(:,j))*(realPOS(:,j)-hatPOS(:,j))';

                my_aPos_state_k(:,j) = hat_v_k_POS(:,j) ;   
                my_aPos_P = (realPOS(:,j)-hat_v_k_POS(:,j))*(realPOS(:,j)-hat_v_k_POS(:,j))';
            end            
            if j >= 4                
                aPri_P = A * aPos_P * A' + Q_k;           
                velocityTerm = hatPOS(:,j) - hatPOS(:,j-1);              
                aPri_state_k(:,j) = aPos_state_k(:,j-1) + (deltaTime * velocityTerm)*10 - mean_w;
            
                % Kalman Gain part !!
                R = v_k * v_k';
                KalmanGain = (aPri_P*H')*pinv(H*aPri_P*H'+R) ;
                
                % Estimation part !!
                aPos_state_k(:,j) = aPri_state_k(:,j) + KalmanGain * (Z_KF-H*aPri_state_k(:,j));
                aPos_P = aPri_P - KalmanGain * H * aPri_P;
                
            %% with Enhanced-KalmanFilter  
                my_aPri_P = A * my_aPos_P * A' + Q_k;               
                my_velocityTerm = hat_v_k_POS(:,j) - hat_v_k_POS(:,j-1);               
                my_aPri_state_k(:,j) = my_aPos_state_k(:,j-1) + (deltaTime * my_velocityTerm)*10 - mean_w;
                
                % Kalman Gain part
                my_R = v_k * v_k';
                my_KalmanGain = (my_aPri_P*H')*pinv(H*my_aPri_P*H'+my_R) ;
                
                % Enhance part
                if s <= 3
                % Estimation part
                    my_aPos_state_k(:,j) = my_aPri_state_k(:,j) + my_KalmanGain * (Z_KF-H*my_aPri_state_k(:,j));
                else
                    my_aPos_state_k(:,j) = my_aPri_state_k(:,j) + my_KalmanGain * (Z-H*my_aPri_state_k(:,j));
                end

                my_aPos_P = my_aPri_P - my_KalmanGain * H * my_aPri_P;

            %% error of filters 
                error_k(:,j) = norm(realPOS(:,j)-hatPOS(:,j));
                error_v_k_k(:,j) = norm(realPOS(:,j)-hat_v_k_POS(:,j));               
                error_kf_k(:,j) = norm(realPOS(:,j)-aPos_state_k(:,j));
                my_error_kf_k(:,j) = norm(realPOS(:,j)-my_aPos_state_k(:,j));
            end
        end
        aPri_state = aPri_state + my_aPri_state_k/iter;
        aPos_state = aPos_state + my_aPos_state_k/iter;
        my_aPri_state = my_aPri_state + my_aPri_state_k/iter;
        my_aPos_state = my_aPos_state + my_aPos_state_k/iter;
        
        error = error + error_k/iter;
        error_v_k = error_v_k + error_v_k_k/iter;
        error_kf = error_kf + error_kf_k/iter;
        my_error_kf = my_error_kf + my_error_kf_k/iter;
    
    end
    ToA.aPri_state_s(2*s-1:2*s,:) = aPri_state;
    ToA.aPos_state_s(2*s-1:2*s,:) = aPos_state;
    ToA.my_aPri_state_s(2*s-1:2*s,:) = my_aPri_state;
    ToA.my_aPos_state_s(2*s-1:2*s,:) = my_aPos_state;
    
    ToA.error_s(s,:) = error;
    ToA.error_v_k_s(s,:) = error_v_k;
    ToA.error_kf_s(s,:) = error_kf;
    ToA.my_error_kf_s(s,:) = my_error_kf;

end

figure;

% graph of error_s
semilogy(log10(sigma_index), mean(ToA.error_s, 2)', 'b'); 
hold on;

% graph of error_v_k_s
semilogy(log10(sigma_index), mean(ToA.error_v_k_s, 2)', 'r-*'); 

% graph of error_kf_s
semilogy(log10(sigma_index), mean(ToA.error_kf_s, 2)', 'c-p'); 

% graph of my_error_kf_s
semilogy(log10(sigma_index), mean(ToA.my_error_kf_s, 2)', 'm-d'); 

% add legend
legend('error_s', 'error_s with v_k', 'error_s with KF', 'error_s with KF');

% add title and labels
title("Performance simulation of ToA using filters");
xlabel('sigma^2');
ylabel('Mean Error');

% add grid
grid on;
