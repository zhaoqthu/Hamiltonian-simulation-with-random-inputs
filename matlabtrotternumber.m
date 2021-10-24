%%Worst-case empirical results are generated here. 
%Compute empirical trotter number
nmax = 5; % largest number of qubis 
num_trials = 5;  %for each n generate 5 instances of Hamiltonians
epsilon = 1e-3;
em_trotter_num = zeros(nmax-2,num_trials+1); % initial 

for n = 3: nmax
    fprintf('n = %d\n', n);
    holder = zeros(1,num_trials);
    init_guess = n * 1e4;

    for counter = 1:num_trials
        [A,B] = heisenberg(n,1); 
        C=0;
        %[A,B,C] = powerlaw(n,4);
        
holder(counter) = empirical_trotter_number(A, B, C, epsilon, init_guess);
em_trotter_num(n-2,counter) = holder(counter);
    end
    
    em_trotter_num(n-2,counter+1) = mean(holder);  % the average value
    fprintf("\n")
end
% % 



%% functions we used
%One-dimension Heisnberg model, h is the random magnetic flied
function [A,B] = heisenberg(n,h)
X = sparse([0 1;1 0]);
Y = sparse([0 -1i;1i 0]);
Z = sparse([1 0;0 -1]);
XX = kron(X,X);
YY = real(kron(Y,Y));
ZZ = kron(Z,Z);

A = sparse(2^n,2^n);
B = sparse(2^n,2^n);

for j = 1:floor(n/2)
    A = A + kron(speye(2^(2*j-2)), kron(XX, speye(2^(n - 2*j)) ) );
    A = A + kron(speye(2^(2*j-2)), kron(YY, speye(2^(n - 2*j)) ) );
    A = A + kron(speye(2^(2*j-2)), kron(ZZ, speye(2^(n - 2*j)) ) );
    A = A + h*(2*rand(1)-1) .* kron(speye(2^(2*j-2)), kron(Z, speye(2^(n - 2*j + 1)) ) );
end

for j = 1:(ceil(n/2)-1)
    B = B + kron(speye(2^(2*j-1)), kron(XX, speye(2^(n - 2*j-1)) ) );
    B = B + kron(speye(2^(2*j-1)), kron(YY, speye(2^(n - 2*j-1)) ) );
    B = B + kron(speye(2^(2*j-1)), kron(ZZ, speye(2^(n - 2*j-1)) ) );
    B = B + h*(2*rand(1)-1) .* kron(speye(2^(2*j-1)), kron(Z, speye(2^(n - 2*j)) ) );
end
end

%One-dimension Heisnberg model with power law interactions, h is the random magnetic flied
% alp=0, 4 are used in our numerical results.
function [A,B,C] = powerlaw(n,alp)
X = sparse([0 1;1 0]);
Y = sparse([0 -1i;1i 0]);
Z = sparse([1 0;0 -1]);
A = sparse(2^n,2^n);
B = sparse(2^n,2^n);
C = sparse(2^n,2^n);

for j= 1:n-1
    h=(2*rand(1)-1);
   for j1=j+1: n
          left = speye(2^(j-1));
          mid=  speye(2^(j1-j-1));
          right =speye(2^(n-j1));
      A= A+ kron(left,kron(X,kron(mid,kron(X,right))))/(j1-j)^alp;
      B= B+ kron(left,kron(Y,kron(mid,kron(Y,right))))/(j1-j)^alp;
      C= C+ kron(left,kron(Z,kron(mid,kron(Z,right))))/(j1-j)^alp;
   end
         left1= speye(2^(j-1));
         right1 =speye(2^(n-j));
      
       C=C+ h.* kron(left1,kron(Z, right1));
 end

end




% worst case empirical PF1 with two terms A, B
function err= empirical_errorW(A,B, r)
%%
[d,~] = size(A);
n = log2(d);
dt = n/r;

U0 = expm(-1i*dt*A) * expm(-1i*dt*B);
U = expm(-1i*n*(A+B));
err = norm(U0^r - U , 2);
end



% worst case empirical PF2 with two terms A, B
function err = empirical_errorW2(A,B, r)
%%
[d,~] = size(A);
n = log2(d);
t=n;
dt = n/r;

U0 = expm(-1i*dt*A/2) * expm(-1i*dt*B)*expm(-1i*dt*A/2) ;
U = expm(-1i*t*(A+B));
err = norm(U0^r - U , 2);
end

% worst case empirical PF1 and PF2 with three terms A, B, C
%ord=1, PF1; ord=2, PF2
function err = empirical_errorWP(A,B,C,r,ord)
%%
[d,~] = size(A);
n = log2(d);
t=n;
dt = n/r;
U2 = expm(-1i*dt*A/2) * expm(-1i*dt*B/2)* expm(-1i*dt*C)*expm(-1i*dt*B/2)*expm(-1i*dt*A/2);
U1=  expm(-1i*dt*A) * expm(-1i*dt*B)* expm(-1i*dt*C);
U0 = expm(-1i*t*(A+B+C));
if ord==1
 err = norm(U1^r - U0 , 2);
else
 err= norm(U2^r - U0 , 2);
end
end

% function of binary search 
function rmin = binary_search(func, epsilon, init_r)
% func is a monotone decreasing function
left = 1;
right = init_r;

while func(right) > epsilon
    fprintf("Initial guess is too small! 10 times the initial guess!");
    right = 10 * init_r;
end

if func(left) <= epsilon
    rmin = 1;
else
    
    while left < right - 1
        mid = ceil((left + right)/2);
        if func(mid) > epsilon
            left = mid;
        else
            right = mid;
        end
     %   fprintf("left:%d, right:%d\n",left,right);
    end
    
    rmin = right;

end

end

% find the minimum number of segments r with different models and methods,
% init_r is the upper bound of initial guess 
function tn = empirical_trotter_number(A, B, C, epsilon, init_r)


 % func_em = @(r) empirical_errorW(A,B,r);  % PF1 
 
  func_em = @(r) empirical_errorW2(A, B, r); %PF2
 %func_em = @(r) empirical_errorWP(A,B,C,r,2); %Power law , PF1 and PF2 
   tn= binary_search(func_em, epsilon, init_r); 
end


