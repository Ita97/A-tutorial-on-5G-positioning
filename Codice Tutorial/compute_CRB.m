function CRB = compute_CRB(u, s, sigma, method)

    H = compute_H(u, s, 2, method);
    CRB = sqrt(sigma^2*trace(eye(2)/(H'*H)));
    
end