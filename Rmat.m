function [R] = Rmat(epps,rho,rrhop,kappa,r1,r2,r3,t,erho1,erho2,erho3,itau, N)
%RMAT Constructs the resitance matrix for body with epps and N
R=zeros(6,6,N+1);

S1=@(s,phi)r1(s)+epps*rho(s).*erho1(s,phi);
S2=@(s,phi)r2(s)+epps*rho(s).*erho2(s,phi);
S3=@(s,phi)r3(s)+epps*rho(s).*erho3(s,phi);

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) [1,0,0],r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,1,:)=F;
R(4:6,1,:)=L;

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) [0,1,0],r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,2,:)=F;
R(4:6,2,:)=L;

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) [0,0,1],r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,3,:)=F;
R(4:6,3,:)=L;

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) cross([1,0,0],[S1(s,phi),S2(s,phi),S3(s,phi)]),r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,4,:)=F;
R(4:6,4,:)=L;

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) cross([0,1,0],[S1(s,phi),S2(s,phi),S3(s,phi)]),r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,5,:)=F;
R(4:6,5,:)=L;

[F,L] = TBT(epps,rho,rrhop,kappa,@(s,phi) cross([0,0,1],[S1(s,phi),S2(s,phi),S3(s,phi)]),r1,r2,r3,t,erho1,erho2,erho3,itau, N);
R(1:3,6,:)=F;
R(4:6,6,:)=L;


end

