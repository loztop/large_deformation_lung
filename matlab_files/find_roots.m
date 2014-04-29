function [z] = find_roots(v)
x = -0.1:0.1:20;

% figure;
%plot(x,besselj(1,x)-((1-v).*x.*besselj(0,x)/(1-2*v))),
%hold on;

%grid on
%f = @(x) tan(x)+5.6*x;
f = @(x) besselj(1,x)-((1-v).*x.*besselj(0,x)/(1-2*v));

N = 30; % Fineness of the search, bigger is finer
z = zeros(1,N);

for jj = 1:20
    for n=1:N
        try
            z(n) = fzero(f,jj+[(n-1)/N n/N]);
            if abs(f(z(n)))>1e-3
                z(n) = 0;
            else
                fprintf('Found a root at: %.5f\n',z(n))
            end
        catch
        end
    end
end

%x = 0:pi/50:10*pi;
y = f(x);

%plot(z,zeros(1,length(z)),'o',x,y,'-')
%line([0 10*pi],[0 0],'color','black')
%axis([0 10*pi -1000.0 3000.0])
z = sort(z(z~=0)); % Return only the needed values 