function [u_k,delta_rho, H] = Non_linear_LS_AoA2(rho, u_0, s, K,stop_cond, dim)
% Perform NLS solution with Iterative search (Gauss-Newton)
% rho: measurements (AOA) in rad, rho(:,1) azimuth, rho(:,2) elevation
% u_0: starting solution
% s: BS positions
% K: number of iterations 
% dim: dimension (2-3 D)

    u_k = u_0(1:dim);
    mu_k = 0.01;
    u_k_prec=u_0+100;
    damp_factor = 0.01;
    for k=1:K
       delta_rho = compute_delta_rho (s, rho, u_k, dim);
       H = compute_H (s, u_k, dim);
       C = (H'*H);
       % ricondiziona matrice mal condizionata
       if rcond(C) <= 1e-4
           %disp("sono entrato");
           C = C + damp_factor*eye(size(C));
       end
       u_k = u_k + (mu_k * C^(-1) * (H') * delta_rho)'; % Gauss-Newton
       error=sqrt(sum((u_k(1:2)-u_k_prec(1:2)).^2));
       
       if error<=stop_cond
           %disp('stop')
           break;
       end
       u_k_prec=u_k;
    end


    function H = compute_H (s, u, dim)
        N = size(s,1);
        d_h = zeros(N,1);
        d = zeros(N,dim); 
        H = zeros(N*(dim-1),dim);
        for i=1:N
            d(i,:) = s(i,1:dim)-u(1:dim); % d(i,1)=d_x, d(i,2)=d_y, d(i,3)=d_z
            d_h(i) = norm(d(i,1:2));
        end
        for i=1:N
            H(i,1) = d(i,2)/(d_h(i))^2; 
            H(i,2) = -d(i,1)/(d_h(i))^2;
        end
        if dim==3
            for i=1:N
                j = i+N;
                d_square = (norm(d(i,:)))^2; % d^2
                H(j,1) = (d(i,3)*d(i,1))/(d_h(i)*d_square); 
                H(j,2) = (d(i,3)*d(i,2))/(d_h(i)*d_square);
                H(j,3) = -d_h(i)/d_square;
            end
        end
    end
    
    function delta_rho = compute_delta_rho (s, rho, u, dim)
        N = size(s,1);
        theta = zeros(N, 1); % azimuth
        phi = zeros(N, 1);   % elevation
        for i=1:N
            theta(i) = atan2((s(i,2)-u(2)),(s(i,1)-u(1)));
        end
        delta_rho = zeros(N*(dim-1), 1);
        for i=1:N
            delta_rho(i) = rho(i,1) - theta(i);
        end
        if dim==3
            for i=1:N
                phi(i) = atan2((s(i,3)-u(3)),((s(i,1)-u(1))*cos(theta(i))+(s(i,2)-u(2))*sin(theta(i))));
                delta_rho(i+N) = rho(i,2) - phi(i);
            end
        end
    end

end 