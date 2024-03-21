function H = compute_H(u,s,dim,method)
   u = u(1:dim);
   s = s(:,1:dim);
   N = size(s,1);
   if strcmpi(method, 'RTT') || strcmpi(method, 'ToA')
        H = zeros(N, dim);
        for i=1:N
            H(i,:) = compute_a(s(i,:), u);
        end
        
   end
   if strcmpi(method, 'TDoA')
        a = zeros(N, dim);
        H = zeros(N-1, dim);
        for i=1:N
            a(i,:) = compute_a(s(i,:), u);
        end
        for i=2:N
            H(i-1,:) = a(1,:) - a(i,:); 
        end
   end
   if strcmpi(method, 'AoA')
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
%         for i=1:N
%             H(i,1) = d(i,2)/d_h(i); 
%             H(i,2) = -d(i,1)/d_h(i);
%         end
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
   if strcmpi(method, 'RTT+AoA')
        d_h = zeros(N,1);
        d = zeros(N, dim); 
        H = zeros(N*dim,dim);
        for i=1:N
            d(i,:) = s(i,1:dim)-u(1:dim); % d(i,1)=d_x, d(i,2)=d_y, d(i,3)=d_z
            d_h(i) = norm(d(i,1:2)); % d_h=sqrt(d_x^2+d_y^2)
        end
        for i=1:N %azimuth
            H(i,1) = d(i,2)/(d_h(i))^2; 
            H(i,2) = -d(i,1)/(d_h(i))^2;
        end
%         for i=1:N %azimuth
%             H(i,1) = d(i,2)/d_h(i); 
%             H(i,2) = -d(i,1)/d_h(i);
%         end
        for i=1:N  % TOA
            H(i+N*(dim-1),:) = compute_a(s(i,:), u); 
        end
        if dim==3  % elevation
            for i=1:N  
                j = i+N;
                d_square = (norm(d(i,:)))^2; % d^2
                H(j,1) = (d(i,3)*d(i,1))/(d_h(i)*d_square); 
                H(j,2) = (d(i,3)*d(i,2))/(d_h(i)*d_square);
                H(j,3) = -d_h(i)/d_square;
            end
        end
   end
   if strcmpi(method, 'TDoA+AoA')
        d_h = zeros(N,1);
        d = zeros(N, dim); 
        a = zeros(N, dim);
        H = zeros(N*dim-1,dim);
        for i=1:N
            d(i,:) = s(i,1:dim)-u(1:dim); % d(i,1)=d_x, d(i,2)=d_y, d(i,3)=d_z
            d_h(i) = norm(d(i,1:2)); % d_h=sqrt(d_x^2+d_y^2)
        end
        for i=1:N
            a(i,:) = (s(i,1:dim)-u(1:dim))./norm(s(i,1:dim)-u(1:dim));
        end
        for i=1:N %azimuth
            H(i,1) = d(i,2)/(d_h(i))^2; 
            H(i,2) = -d(i,1)/(d_h(i))^2;
        end
%         for i=1:N %azimuth
%             H(i,1) = d(i,2)/d_h(i); 
%             H(i,2) = -d(i,1)/d_h(i);
%         end
        for i=2:N  % TDOA
            H(i-1+N*(dim-1),:) = a(1,:) - a(i,:); 
        end
        if dim==3  % elevation
            for i=1:N  
                j = i+N;
                d_square = (norm(d(i,:)))^2; % d^2
                H(j,1) = (d(i,3)*d(i,1))/(d_h(i)*d_square); 
                H(j,2) = (d(i,3)*d(i,2))/(d_h(i)*d_square);
                H(j,3) = -d_h(i)/d_square;
            end
        end
   end
    
   function a = compute_a(s, u)
         d_i = norm(s-u);
         a = (s-u)./d_i;
   end
end