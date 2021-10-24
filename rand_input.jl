
module rand_input
export RandInp,triangle_errorpower2, search_trotter_num_triangle23, trotteran43, trotteran63, trotteran23, init_heisenberg!,trotteran4,trotteran6,search_trotter_num_triangle3, triangle_errorpcount,triangle_errorpower, triangle_errorSe,init_powerlaw!, trotteranp, search_empir_numR,trotteran,trotteran2, trotter_state, sp_trace, binary_search, interference_error, triangle_error, search_trotter_num_inteference, search_trotter_num_triangle

using LinearAlgebra
using SparseArrays
using ExpmV
using QuantumInformation
using Statistics
using Cubature
using Arpack




# generate different types of random input
#Style: 1 local random; 2 random computational basis; other, haar random
# n number of qubit,
function RandInp(n, style)
 X=[0 1;1 0]
 d=2^n;
if style==1
d1=CUE(2)
Inp=1;
for j in 1:n
  U1=rand(d1);
  psi=U1*ket(1,2)
  Inp=kron(Inp,psi);
end
elseif style==2
  Inp=1;
  for j in 1:n
    a=rand(0:1);
    psi=[1; 0]
    psi=X^a*psi;
    Inp=kron(Inp,psi);
   end
else
   se=HaarKet{2}(d);
   Inp=rand(se);
end
return Inp
end

# Heisenberg model Hamiltonian
# s is the random magnetic field
function init_heisenberg!(A::AbstractSparseArray, B::AbstractSparseArray, n::Int64,s)

    # Preparing building blocks
    XX = sparse([0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]);
    YY = sparse([0 0 0 -1;0 0 1 0;0 1 0 0;-1 0 0 0]);

    Z = sparse([1 0;0 -1]);
    ZZ = kron(Z,Z);
    Ide= sparse([1 0;0 1]);

    # Building A
    for j in 1:Int(floor(n/2))
        rh=2*rand(1)-[1];
        h=rh*s;
        if j == 1
            left = 1;
            right = sparse(Matrix(1*I, 2^(n-2), 2^(n-2)));
        elseif n == 2*j
            left = sparse(Matrix(1*I, 2^(n-2), 2^(n-2)));
            right = 1;
        else
            left = sparse(Matrix(1*I, 2^(2*j-2), 2^(2*j-2)));
            right = sparse(Matrix(1*I, 2^(n-2*j), 2^(n-2*j)));
        end

        A .+= kron(left, XX, right);
        A .+= kron(left, YY, right);
        A .+= kron(left, ZZ, right);
        A .+= h.*kron(left, Z, Ide, right);

    end

    # Building B
    for j in 1:Int(ceil(n/2 - 1))
        rh=2*rand(1)-[1];
        h=rh*s;
        if n == 2*j+1
            left = sparse(Matrix(1*I, 2^(n-2), 2^(n-2)));
            right = 1;
        else
            left = sparse(Matrix(1*I, 2^(2*j-1), 2^(2*j-1)));
            right = sparse(Matrix(1*I, 2^(n-2*j-1), 2^(n-2*j-1)));
        end

        B .+= kron(left, XX, right);
        B .+= kron(left,YY,right);
        B .+= kron(left,ZZ,right);
       B .+= h.*kron(left, Z, Ide, right);

    end
end

# Heisenberg model Hamiltonian with power law interaction
# s is the random magnetic field
function init_powerlaw!(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray, n::Int64, alp)

    # Preparing building blocks
    XX = sparse([0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]);
    YY = sparse([0 0 0 -1;0 0 1 0;0 1 0 0;-1 0 0 0]);
    Z = sparse([1 0;0 -1]);
    ZZ = kron(Z,Z);
    X= sparse([0 1;1 0]);
    Y= sparse([0 -1im;1im 0]);
    count=0;
    # Building A B C
    for j in 1:n-1
        h=2*rand(1)-[1];
        for j1 in j+1: n
            left = sparse(Matrix(1*I, 2^(j-1), 2^(j-1)));;
            mid=  sparse(Matrix(1*I, 2^(j1-j-1), 2^(j1-j-1)));
            right = sparse(Matrix(1*I, 2^(n-j1), 2^(n-j1)));

        A .+= kron(left, X, mid, X, right)./(j1-j)^alp;
        B .+= kron(left, Y, mid, Y, right)./(j1-j)^alp;
        C .+= kron(left, Z, mid, Z, right)./(j1-j)^alp;
         end

        left1 = sparse(Matrix(1*I, 2^(j-1), 2^(j-1)));;
        right1 = sparse(Matrix(1*I, 2^(n-j), 2^(n-j)));

        C.+=h.* kron(left1, Z, right1);

     end
end


#1.  Generate empirical results

# Average error with a fixed number of segments r, PF1 ,  H=A+B
# samp: sample number for random inputs
# style: type of random inputs
# type: 1. average trotter error, 2. standard deviation
function trotteran(A::AbstractSparseArray, B::AbstractSparseArray, t, r, samp, style, type)
dt = t/r;
H= A + B;
d, ~ = size(H);
n=Int(log2(d));
sample=samp;
ep=zeros(sample);

for j in 1:sample
psi1=RandInp(n, style);
a2=expmv(-1im .* t, H, psi1);
storage = zeros(Complex{Float64}, d);
storage.=psi1;
for j in 1:r
   storage .= expmv(-1im .* dt, A, storage);
   storage .= expmv(-1im .* dt, B, storage);
end
   a4=storage;
 ep[j]=norm(a4-a2,2);
end
v=stdm(ep,mean(ep))
if type==1
return  mean(ep)
  else
return v
end
end

# Average error with a fixed number of segments r,  PF1 ,  H=A+B+C
function trotteranp(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray, t, r, samp,type)
dt = t/r;
H= A + B + C;
d, ~ = size(H);
sample=samp;
ep=zeros(sample);
se=HaarKet{2}(d);
for j in 1:sample
psi1=rand(se);
a2=expmv(-1im .* t, H, psi1);
storage = zeros(Complex{Float64}, d);
storage.=psi1;
for j in 1:r
   storage .= expmv(-1im .* dt, A, storage);
   storage .= expmv(-1im .* dt, B, storage);
   storage .= expmv(-1im .* dt, C, storage);
end
   a4=storage;
 ep[j]=norm(a4-a2,2);
end
v=stdm(ep,mean(ep))
if type==1
return  mean(ep)
  else
return v
end
end

# Average error with a fixed number of segments r, PF2 ,  H=A+B
function trotteran2(A::AbstractSparseArray, B::AbstractSparseArray, t, r, samp,type)
    dt = t/r;
    H= A + B;
    d, ~ = size(H);
    sample=samp;
    ep=zeros(sample);
    va=0;
    se=HaarKet{2}(d);
    for j in 1:sample
    psi1=rand(se);
    a2=expmv(-1im .* t, H, psi1);
    storage = zeros(Complex{Float64}, d);
    storage.=psi1;
    for j in 1:r
     storage .=trotter_state(A,B,dt,storage);
    end
       a4=storage;
     ep[j]=norm(a4-a2,2);
    end
    v=stdm(ep,mean(ep))
    if type==1
    return  mean(ep)
      else
    return v
    end
    end

# Average error with a fixed number of segments r,  PF2 ,  H=A+B+C
function trotteran23(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray,t, r, samp)
                dt = t/r;
                H= A + B+ C;
                d, ~ = size(H);
                sample=samp;
                ep=zeros(sample);
                se=HaarKet{2}(d);
                for j in 1:sample
                psi1=rand(se);
                a2=expmv(-1im .* t, H, psi1);
                storage = zeros(Complex{Float64}, d);
                storage.=psi1;
                for j in 1:r
                storage .=trotter_state23(A,B,C,dt,storage);
                end
                   a4=storage;
                 ep[j]=norm(a4-a2,2);
                end

                return  mean(ep)
    end


# Average error with a fixed number of segments r,  PF4 ,  H=A+B
function trotteran4(A::AbstractSparseArray, B::AbstractSparseArray, t, r, samp)
        dt = t/r;
        H= A + B;
        d, ~ = size(H);
        sample=samp;
        ep=zeros(sample);
        se=HaarKet{2}(d);
        for j in 1:sample
        psi1=rand(se);
        a2=expmv(-1im .* t, H, psi1);
        storage = zeros(Complex{Float64}, d);
        storage.=psi1;
        p=1/(4-4^(1/3));
        for j in 1:r
         storage .=trotter_state4(A,B,dt,storage);
        end
           a4=storage;
         ep[j]=norm(a4-a2,2);
        end

        return  mean(ep)
        end


# Average error with a fixed number of segments r,  PF4 ,  H=A+B+C
function trotteran43(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray,t, r, samp)
        dt = t/r;
        H= A + B+C;
        d, ~ = size(H);
        sample=samp;
        ep=zeros(sample);
        se=HaarKet{2}(d);
        for j in 1:sample
          psi1=rand(se);
        a2=expmv(-1im .* t, H, psi1);
        storage = zeros(Complex{Float64}, d);
        storage.=psi1;
        for j in 1:r
            storage .=trotter_state43(A,B,C,dt,storage);
        end
        a4=storage;
        ep[j]=norm(a4-a2,2);
        end
        return  mean(ep)
        end



# Average error with a fixed number of segments r,  PF6 ,  H=A+B
function trotteran6(A::AbstractSparseArray, B::AbstractSparseArray, t, r, samp)
            dt = t/r;
            H= A + B;
            d, ~ = size(H);
            sample=samp;
            ep=zeros(sample);
            se=HaarKet{2}(d);
            for j in 1:sample
            psi1=rand(se);
            a2=expmv(-1im .* t, H, psi1);
            storage = zeros(Complex{Float64}, d);
            storage.=psi1;
            p=1/(4-4^(1/5));
            for j in 1:r
            storage .=trotter_state6(A,B,dt,storage);
            end
               a4=storage;
             ep[j]=norm(a4-a2,2);
            end

            return  mean(ep)
end

# Average error with a fixed number of segments r,  PF6 ,  H=A+B+C
function trotteran63(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray,t, r, samp)
            dt = t/r;
            H= A + B+C;
            d, ~ = size(H);
            sample=samp;
            ep=zeros(sample);
            se=HaarKet{2}(d);
            for j in 1:sample
            psi1=rand(se);
            a2=expmv(-1im .* t, H, psi1);
            storage = zeros(Complex{Float64}, d);
            storage.=psi1;
            for j in 1:r
            storage .=trotter_state63(A,B,C,dt,storage);
            end
               a4=storage;
             ep[j]=norm(a4-a2,2);
            end

            return  mean(ep)
end



# output the state after t time evolution with PF2, and H=A+B
function trotter_state(A::AbstractSparseArray, B::AbstractSparseArray, t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
    #dt = t/r;
#    println("storage")
        storage .= expmv(-1im .* t, A/2, storage);
        storage .= expmv(-1im .* t, B, storage);
        storage .= expmv(-1im .* t, A/2, storage);
    return storage
end


# output the state after t time evolution with PF2, and H=A+B+C
function trotter_state23(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray,t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
    #dt = t/r;
#    println("storage")
        storage .= expmv(-1im .* t, A/2, storage);
        storage .= expmv(-1im .* t, B/2, storage);
        storage .= expmv(-1im .* t, C, storage);
        storage .= expmv(-1im .* t, B/2, storage);
        storage .= expmv(-1im .* t, A/2, storage);
    return storage
end




function trotter_state4(A::AbstractSparseArray, B::AbstractSparseArray, t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
p=1/(4-4^(1/3));
    storage .=trotter_state(A,B,p*t,storage);
    storage .=trotter_state(A,B,p*t,storage);
    storage .=trotter_state(A,B,(1-4*p)*t,storage);
    storage .=trotter_state(A,B,p*t,storage);
    storage .=trotter_state(A,B,p*t,storage);
    return storage
end

# output the state after t time evolution with PF4, and H=A+B+C
function trotter_state43(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray,t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
p=1/(4-4^(1/3));
    storage .=trotter_state23(A,B,C,p*t,storage);
    storage .=trotter_state23(A,B,C,p*t,storage);
    storage .=trotter_state23(A,B,C,(1-4*p)*t,storage);
    storage .=trotter_state23(A,B,C,p*t,storage);
    storage .=trotter_state23(A,B,C,p*t,storage);
    return storage
end

# output the state after t time evolution with PF6, and H=A+B
function trotter_state6(A::AbstractSparseArray, B::AbstractSparseArray, t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
p=1/(4-4^(1/5));
    storage .=trotter_state4(A,B,p*t,storage);
    storage .=trotter_state4(A,B,p*t,storage);
    storage .=trotter_state4(A,B,(1-4*p)*t,storage);
    storage .=trotter_state4(A,B,p*t,storage);
    storage .=trotter_state4(A,B,p*t,storage);
    return storage
end



# output the state after t time evolution with PF6, and H=A+B+C
function trotter_state63(A::AbstractSparseArray, B::AbstractSparseArray, C::AbstractSparseArray, t, psi0)
   d, ~ = size(A);
   storage = zeros(Complex{Float64}, d);
    storage .= psi0;
    p=1/(4-4^(1/5));
    storage .=trotter_state43(A,B,C,p*t,storage);
    storage .=trotter_state43(A,B,C,p*t,storage);
    storage .=trotter_state43(A,B,C,(1-4*p)*t,storage);
    storage .=trotter_state43(A,B,C,p*t,storage);
    storage .=trotter_state43(A,B,C,p*t,storage);
    return storage
end

#binary search function with range [init_left, init_right],
#approximated search, the search ends when left > right- step
function binary_search(func, epsilon, init_left, init_right, step)
    left = init_left;
    right = init_right;

    while func(right) > epsilon
        println("Initial guess is too small! 10 times the initial guess!");
        right *= 10;
    end

    if func(left) <= epsilon
        rmin = left;
    else
        while left < right - step
            mid = ceil((left + right)/2);
            if func(mid) > epsilon
                left = mid;
            else
                right = mid;
            end
            println("left: ", left, " right: ", right);
        end

        rmin = right;
    end

    return rmin
end

#find the minimum r
function search_empir_numR(A, B, C,t, ϵ, iniguessL, iniguessR, step)
    H = A + B + C;
    d, ~ = size(H);
#change it to other random error function, for PF1, 2,4,6 ,power law
#empi_error(r) =trotteranp(A,B,C, t,r,20,1);

#empi_error(r) =trotteran23(A,B,C, t,r,20);
#empi_error(r) =trotteran43(A,B,C, t,r,20);
#empi_error(r) =trotteran6(A,B,t,r,20);
#empi_error(r) =trotteran4(A,B,t,r,20);
#empi_error(r) =trotteran2(A,B,t,r,20,1);
empi_error(r) =trotteran(A,B,t,r,20,0,1);
    init_left = iniguessL;
    init_right = iniguessR;
    rmin = binary_search(empi_error, ϵ, init_left, init_right, step)
    #va=trotteranp(A,B,C,t,rmin,20,2);
    va=trotteran(A,B,t,rmin,20,0,2);
 # corresponding standard deviation with rmin segments
 #change it to other random error function,
  return  [rmin,va]
end












# calculate trace for sparse matrix
function sp_trace(A::AbstractSparseArray)
    m,n = size(A);
    if m != n
        println("Input must be a square matrix. Operation Aborted.");
        return nothing
    else
        tr = 0;
        for i in 1:n
            tr += A[i,i]
        end
        return tr
    end
end

# output commutator
function commutator!(storage::AbstractSparseArray, A::AbstractSparseArray, B::AbstractSparseArray)
    storage .= A * B - B * A;
end


# Interference error bound

function interference_error(A::AbstractSparseArray, B::AbstractSparseArray, t, r)
    H = A + B;
    d, ~ = size(H);
    n = log2(d);
    dt = t/r;
    C = 2048/(ℯ^2 * (ℯ-1));
    z=svds(H,nsv = 1)[1];
    Hor=z.S[1];
    tHor = Hor*dt;

    nest_commutator = spzeros(d, d);
    commutator!(nest_commutator, H, B);
    commutator!(nest_commutator, A, nest_commutator)
    com_AB = spzeros(d, d);
    commutator!(com_AB, A, B)
    zab=svds(com_AB,nsv = 1)[1];
    sp_AB=zab.S[1];
    term_a = 1/(1 - ℯ*tHor)^2 - 1;
    term_b = sp_trace(nest_commutator * nest_commutator') * dt^6/36;
    trVV = d * C^2 * ℯ^2 * Hor^2 * dt^6 * term_a/(36 * n) + term_b;
    Σ = sqrt(d*n) * dt^2/(1 - tHor)^2;
    sqrt_trΔ = 2*r*Σ/t + r * sqrt(trVV);
    M1_term = 1/(1 - r * 0.5*dt^2*sp_AB)
    if r * 0.5*dt^2*sp_AB>1
        ep=1;
    else ep=sqrt_trΔ * M1_term/sqrt(d)
    end
    return ep
end


function search_trotter_num_inteference(A, B, t, ϵ)
    intfr_error(r) = interference_error(A, B, t, r);
    H = A + B;
    init_left = 1*10^3;
    init_right = 3*10^5;
    rmin = binary_search(intfr_error, ϵ, init_left, init_right,1)
end







# triangle bound

# triangle error bound PF1
function triangle_error(A::AbstractSparseArray, B::AbstractSparseArray, t, r)
    H = A + B;
    d, ~ = size(H);
    n = log2(d);
    dt = t/r;

    com_AB = spzeros(d, d);
    commutator!(com_AB, A, B)

    term1 = sp_trace(com_AB * com_AB');
    step_err = 0.5 * dt^2 * sqrt(term1/d);

    return r * step_err
end

# triangle error bound PF2
function triangle_errorSe(A::AbstractSparseArray, B::AbstractSparseArray, t, r)
    H = A + B;
    d, ~ = size(H);
    n = log2(d);
    dt = t/r;

    com_AB = spzeros(d, d);
    com_ABB= spzeros(d, d);
    com_AAB= spzeros(d, d);
    commutator!(com_AB, A, B)
    commutator!(com_ABB, com_AB, B)
    commutator!(com_AAB, com_AB, A)

    term1 = sp_trace(com_ABB * com_ABB');
    term2 = sp_trace(com_AAB * com_AAB');
    step_err = dt^3 * sqrt(term1/d)/12+ dt^3 * sqrt(term2/d)/24;

    return r * step_err
end

function search_trotter_num_triangle(A, B, t, ϵ)
    trig_error(r) = triangle_error(A, B, t, r);  #PF1
    #trig_error(r) = triangle_errorSe(A, B, t, r); #PF2

    init_left = 1;
    init_right = 1e3;
    rmin = binary_search(trig_error, ϵ, init_left, init_right,1 )
end


# H=A+B+C case

# PF1
function triangle_errorpower(A::AbstractSparseArray, B::AbstractSparseArray,C::AbstractSparseArray)
    H = A + B +C;
    d, ~ = size(H);
    n = log2(d);

    com_BC = spzeros(d, d);
    com_ABC= spzeros(d, d);

    commutator!(com_BC, B, C)
    commutator!(com_ABC, A, B+C)

    term1 = sp_trace(com_ABC * com_ABC');
    term2 =sp_trace(com_BC * com_BC');
  err= sqrt(term1/d)+sqrt(term2/d)
    return err
end

# PF2
function triangle_errorpower2(A::AbstractSparseArray, B::AbstractSparseArray,C::AbstractSparseArray)
    H = A + B +C;
    d, ~ = size(H);
    n = log2(d);
    com_BC = spzeros(d, d);
    com_AABC= spzeros(d, d);
    com_ABCBC= spzeros(d, d);
    com_BBC= spzeros(d, d);
    com_BCC= spzeros(d, d);
    com_ABC= spzeros(d, d);
    #com_AAB= spzeros(d, d);
    commutator!(com_BC, B, C)
    commutator!(com_BBC, com_BC, B)
    commutator!(com_BCC, com_BC, C)
    commutator!(com_ABC, A, B+C)
    commutator!(com_AABC, A, com_ABC)
    commutator!(com_ABCBC, B+C, com_ABC)

    term1 = sp_trace(com_BCC * com_BCC');
    term2 =sp_trace(com_ABCBC * com_ABCBC');
    term3 =sp_trace(com_BBC * com_BBC');
    term4 =sp_trace(com_AABC * com_AABC');
 err= sqrt(term1/d)/12+sqrt(term2/d)/12+sqrt(term3/d)/24+sqrt(term4/d)/24;
    return err
end

# r PF1 H=A+B+C
function search_trotter_num_triangle3(A, B, C, t, ϵ)
  rmin = triangle_errorpower(A, B, C)*t^2/(2*ϵ)
end


# r PF2 H=A+B+C
function search_trotter_num_triangle23(A, B, C, t, ϵ)
  rmin = triangle_errorpower2(A, B, C)*t^3/ϵ;
  rmin= sqrt(rmin);
    return rmin
end





# counting bound for Powerlaw
function triangle_errorpcount(t,n,tol, alp)
T=0;
TC=0;
TN=0;
truc=1000;
for  i in 1: n-1
    for  j in i+1: n
# here 38400 is truncated, otherwise, it will exceed the max length
        d1=min((i-j)^alp,38400)
        TC=TC+ 8/(d1^4);
        TN=TN+ 8/(d1^2);
        #i1=i
         for j1 in i+1:n
            d2=min((j1-i)^alp,38400)
            T=T+4/(d1^2*d2^2);
         end
         for j1 in j+1: n
            d2=min((j1-j)^alp,38400)
            T=T+4/(d1*d2)^2;
         end
         if i==1
             T=T+0;
         else
            for i1 in 1:i-1
            d2=min((i-i1)^alp,38400)
            T=T+4/(d1*d2)^2;
            end
        end
         for i1 in 1:j-1
            d2=min((i1-j)^alp,38400)
            T=T+4/(d1*d2)^2;
         end
    end
end
  return  (sqrt(2)+1)*sqrt((T-TC)+TN)*0.5*t^2/tol
end







end
