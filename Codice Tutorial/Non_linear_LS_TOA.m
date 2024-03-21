% TOA nls
function [u_k, delta_rho,H] = Non_linear_LS_TOA(rho, u_0, s, K,stop_cond, dim)
% Perform NLS solution with Iterative search (Gauss-Newton)

    u_k = u_0(1:dim);
    s = s(:,1:dim);
    mu_k = 0.01;
    damp_factor = 0.01;
    u_k_prec=u_0+100;
    s = s(:,1:dim);
    for k=1:K
       delta_rho = compute_delta_rho (s, rho, u_k);
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
        a = zeros(N, dim);
        H = zeros(N, dim);
        for i=1:N
            H(i,:) = compute_a(s(i,:), u);
        end
    end

    function a = compute_a(s, u)
         d = norm(s-u);
         a = (u-s)./d;
    end
    
    function delta_rho = compute_delta_rho (s, rho, u)
        N = size(s,1);
        d = zeros(1, N);
        for i=1:N
            d(i) = norm(s(i,:)-u);
        end
        delta_rho = zeros(N, 1);
        for i=1:N
            delta_rho(i) = rho(i) - d(i);
        end
    end

end 