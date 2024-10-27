 function z  = calculate_channel_polarization( epsilon, N)
           z = zeros(N,1);
           z(1) = 1/exp(epsilon); 
           for j =1:log2(N)
               u=2^j;
               count = 1;
               for t=0:u/2-1
                   T =z(count);
                   z(count) = 2*T-T^2;
                   z(u/2+count) = T^2;
                   count=count+1;
               end
               
           end
 end