using LinearAlgebra
using SparseArrays
using ExpmV
using Plots
using QuantumInformation
using Statistics
using Arpack

include("rand_input.jl")

using Main.rand_input


#t=1;
#r=100;
#n=10
#dim = 2^n;
#d = dim;
#A = spzeros(Int64, d, d);
#B = spzeros(Int64, d, d);
#init_heisenberg!(A, B, n);
#er=trotteran(A,B,t,r,30);
tol=1e-3;
ma=3;    #  maximum number of qubits
rm=5;    # sample random magnetic field
#r=zeros(ma);
r0=zeros(ma-2,rm);
r1=zeros(ma-2,rm);
r2=zeros(ma-2,rm);
r3=zeros(ma-2,rm);




for j in 1: ma-2
n=2+j;
println(n)
  for j1 in 1:rm
   t=n;
   dim = 2^n;
   d = dim;
 A = spzeros(Float64, d, d);
 B = spzeros(Float64, d, d);
 C = spzeros(Float64, d, d);
   init_heisenberg!(A, B, n,1);
   #init_powerlaw!(A,B,C,n,0)

# empirical
@time r0[j,j1]=search_empir_numR(A,B,C,t,tol,1e+3,2e+4,1)[1];
#change the function in search_empir_numR


# triangle bound
#r1[j]=search_trotter_num_triangle(A,B,t,tol)
# search_trotter_num_triangle3  for A, B, C, PF1
#search_trotter_num_triangle23 for A, B, C, PF2


# interference bound
# r2[j]= search_trotter_num_inteference(A,B,t,tol)

#Couting bound for power law
#r3[j]= triangle_errorpcount(t,n,tol,4)

end
#mean(r1[j,:])
println(r3)
#println(ev2)
end
#a=zeros(ma);
#b=zeros(ma);


#for j in 1:ma

#a[j] =floor(Int64, mean(r2[j,:]));
#b[j] =stdm(r2[j,:],a[j])
#end


#@time epp[j]=trotteran(A,B,t,r,20);

#trotteran4
#epi[j]= triangle_error(A,B,t,r);
#println(ep[j])
#7e+4,9e+4
#@time r[j]=search_empir_numR(A,B,tem,tol,100,200);
#r2[j]=search_trotter_num_triangle3(A,B,C,tem,tol);
#println(r[j])
#println(epp[j])
#end


#H=A+B;
#r1=search_trotter_num_triangle(A, B, t, tol);
#r2=search_trotter_num_inteference(A, B, t, tol);

#@time   =trottervec(A,B,t,)
#@time a4=trotteran(A, B, t, 200, 20)
#@time a5=trotteran(A, B, t, 27400, 20)
#println(a5)
#@time r3=search_empir_numR(A,B,t,tol,1700,2400);
#@time r3=search_empir_numR(A,B,t,tol,100,0.5e+4);
