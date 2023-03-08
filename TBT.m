function [F,L] = TBT(epps,rho,rrhop,kappa,V,r1,r2,r3,t,erho1,erho2,erho3,itau, N)
%ACSBT Calculates the results of TBT for a prolate ellipsoid
%  INPUTS
%  epps     -   slenderness parameter
%  rho(s)   -   cross sectional shape
%  rrhop(s) -   derivative of cross sectional shape times (rho(s)rho'(s))
%  kappa(s) -   curvature of centerline
%  V(s,phi) -   Velocity on the surface
%  r1(s)    -   centerline in x direction
%  r2(s)    -   centerline in y direction
%  r3(s)    -   centerline in z direction
%  t(s)     -   tangent to centerline
%  erho1(s,phi)-   local radial vector component in x
%  erho1(s,phi)-   local radial vector component in y
%  erho1(s,phi)-   local radial vector component in z
%  itau(s)  -   integrated torsion along the body
%  N        -   Order to expand the solution to
% OUTPUTS
%  F        -   Force from the surface velocity
%  L        -   Torque from the surface velocity


%surface

S1=@(s,phi)r1(s)+epps*rho(s).*erho1(s,phi);
S2=@(s,phi)r2(s)+epps*rho(s).*erho2(s,phi);
S3=@(s,phi)r3(s)+epps*rho(s).*erho3(s,phi);

%prep terms
F= zeros(3,N+1);
L= zeros(3,N+1);

steps = 30;
ds = 2/(steps); %distance between consecutive s points.
s = linspace(-1+(ds/2),1-(ds/2),steps); %s points to be used for colocation method
dphi = 2*pi/(steps+2);
phi=linspace(-pi+dphi/2,pi-dphi/2,steps+2);

[S,Phi]=meshgrid(s,phi);

% geometery and velocity
u=zeros(steps+2,steps);
c=zeros(steps+2,steps);
a=zeros(steps+2,steps);
alpha=zeros(steps+2,steps);
V1=zeros(steps+2,steps);
V2=zeros(steps+2,steps);
V3=zeros(steps+2,steps);
tv=zeros(steps+2,steps,3);
for i=1:steps
    for j=1:(steps+2)
        temp=V(s(i),phi(j));
        V1(j,i)=temp(1);
        V2(j,i)=temp(2);
        V3(j,i)=temp(3);
        tv(j,i,:)=t(s(i));
        rt= rho(s(i))^2 +sqrt(rho(s(i))^4+4*rrhop(s(i))^2);
        c(j,i)=sqrt(rt/2);
        u(j,i) =-2*rrhop(s(i))/rt;
        a(j,i) = 1-epps*rho(s(i))*kappa(s(i))*cos(phi(j)-itau(s(i)));
        %a(j,i) = 1-epps*rho(s(i))*kappa(s(i))*cos(phi(j));
        alpha(j,i) = epps*c(j,i)/a(j,i);
    end
end


%construct zeta_1' and zeta_2' at points s
z1p = Z1p(a,alpha,u);
z2p = Z2p(a,alpha,u);

%construct <1/zeta_1'>  <1/zeta_2'> at points s
z1pib= dphi*sum(1./z1p);
z2pib= dphi*sum(1./z2p);

% Construction of matrix B %slowest part of process
Bl=zeros(3*steps);
Bn=zeros(3*steps);
for i=1:steps
    Bl(((1:3)+(i-1)*3),((1:3)+(i-1)*3))=t(s(i)).*t(s(i))'/z1pib(i) +(eye(3)-t(s(i)).*t(s(i))')/z2pib(i);
    for j=1:(steps)
        Bi =integral(@(s2) 1./sqrt((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)),s(j)-ds/2,s(j)+ds/2);
        B11 =integral(@(s2) ((r1(s(i))-r1(s2)).^2)./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        B12 =integral(@(s2) ((r1(s(i))-r1(s2)).*(r2(s(i))-r2(s2)))./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        B13 =integral(@(s2) ((r1(s(i))-r1(s2)).*(r3(s(i))-r3(s2)))./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        B22 =integral(@(s2) ((r2(s(i))-r2(s2)).*(r2(s(i))-r2(s2)))./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        B23 =integral(@(s2) ((r2(s(i))-r2(s2)).*(r3(s(i))-r3(s2)))./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        B33 =integral(@(s2) ((r3(s(i))-r3(s2)).*(r3(s(i))-r3(s2)))./(((r1(s(i))-r1(s2)).^2+(r2(s(i))-r2(s2)).^2+(r3(s(i))-r3(s2)).^2 +(epps^2).*(rho(s(i)).^2 +rho(s2).^2)).^(3/2)),s(j)-ds/2,s(j)+ds/2);
        Bn(((1:3)+(i-1)*3),((1:3)+(j-1)*3))=Bi*eye(3) +[B11 B12 B13; B12 B22 B23; B13 B23 B33];
    end
end
B=Bl+Bn;


% zeroth itertion
Mival=zeros(3,steps+2,steps);
tV=tv(:,:,1).*V1+tv(:,:,2).*V2+tv(:,:,3).*V3;
Mival(1,:,:)=8*pi*(V1./z2p +tv(:,:,1).*tV.*(1./z1p -1./z2p));
Mival(2,:,:)=8*pi*(V2./z2p +tv(:,:,2).*tV.*(1./z1p -1./z2p));
Mival(3,:,:)=8*pi*(V3./z2p +tv(:,:,3).*tV.*(1./z1p -1./z2p));

Q= reshape(dphi*sum(Mival,2),3*steps,1);

A0=Bl*Q;

fb=B\A0;

MT=reshape(A0- Bl*fb,3,steps);

f0 = zeros(3,steps+2,steps);
l0 = zeros(3,steps+2,steps);
for i=1:(steps+2)
    for j=1:(steps)
        f0(:,i,j)= Mival(:,i,j) -(t(s(j)).*t(s(j))')*MT(:,j)/z1p(i,j) -((eye(3)-t(s(j)).*t(s(j))'))*MT(:,j)/z2p(i,j);
        l0(:,i,j)= cross([S1(S(i,j),Phi(i,j)),S2(S(i,j),Phi(i,j)),S3(S(i,j),Phi(i,j))],f0(:,i,j));
    end
end

f= f0;

F(1,1)= sum(f0(1,:,:),'all')*ds*dphi;
F(2,1)= sum(f0(2,:,:),'all')*ds*dphi;
F(3,1)= sum(f0(3,:,:),'all')*ds*dphi;
L(1,1)= sum(l0(1,:,:),'all')*ds*dphi;
L(2,1)= sum(l0(2,:,:),'all')*ds*dphi;
L(3,1)= sum(l0(3,:,:),'all')*ds*dphi;


%consider replacing this with a while statement

if N~=0
    
    fxi =@(y,x) interp2(phi,s,reshape(f0(1,:,:),steps+2,steps)',x,y,'spline');
    fyi =@(y,x) interp2(phi,s,reshape(f0(2,:,:),steps+2,steps)',x,y,'spline');
    fzi =@(y,x) interp2(phi,s,reshape(f0(3,:,:),steps+2,steps)',x,y,'spline');
    q1m=8*pi*V1;
    q2m=8*pi*V2;
    q3m=8*pi*V3;
    
    %Integration of the force
    
    
    
    for n=1:N
        
        B1v=zeros((steps+2)*steps,1);
        B2v=zeros((steps+2)*steps,1);
        B3v=zeros((steps+2)*steps,1);
        %parfor (l=0:((steps+2)*steps-1),3)
            
        parfor (l=0:((steps+2)*steps-1) ,3)
            %tic
            j=floor(l/(steps+2))+1;
            i=rem(l,steps+2)+1;
            s = linspace(-1+(ds/2),1-(ds/2),steps);
            phi=linspace(-pi+dphi/2,pi-dphi/2,steps+2);
            pt=phi(i);
            st=s(j);
            Ci=@(s2,phi2)1./sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2);
            C11=@(s2,phi2) (S1(st,pt)-S1(s2,phi2)).*(S1(st,pt)-S1(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            C12=@(s2,phi2) (S1(st,pt)-S1(s2,phi2)).*(S2(st,pt)-S2(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            C13=@(s2,phi2) (S1(st,pt)-S1(s2,phi2)).*(S3(st,pt)-S3(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            C22=@(s2,phi2) (S2(st,pt)-S2(s2,phi2)).*(S2(st,pt)-S2(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            C23=@(s2,phi2) (S2(st,pt)-S2(s2,phi2)).*(S3(st,pt)-S3(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            C33=@(s2,phi2) (S3(st,pt)-S3(s2,phi2)).*(S3(st,pt)-S3(s2,phi2))./(sqrt((S1(st,pt)-S1(s2,phi2)).^2+(S2(st,pt)-S2(s2,phi2)).^2+(S3(st,pt)-S3(s2,phi2)).^2).^3);
            BI1 =@(x,y) (Ci(x,y)+C11(x,y)).*fxi(x,y) + C12(x,y).*fyi(x,y) + C13(x,y).*fzi(x,y);
            BI2 =@(x,y) (C12(x,y)).*fxi(x,y) + (Ci(x,y)+C22(x,y)).*fyi(x,y) + C23(x,y).*fzi(x,y);
            BI3 =@(x,y) (C13(x,y)).*fxi(x,y) + C23(x,y).*fyi(x,y) + (Ci(x,y)+C33(x,y)).*fzi(x,y);
            
            a2=a(i,j);
            t12= tv(i,j,1);
            t22= tv(i,j,2);
            t32= tv(i,j,3);
            r12=r1(st);
            r22=r2(st);
            r32=r3(st);
            c2=c(i,j);
            u2=u(i,j);
            r1e=@(s2) a2.*(s2-u2).*t12 +r12;
            r2e=@(s2) a2.*(s2-u2).*t22 +r22;
            r3e=@(s2) a2.*(s2-u2).*t32 +r32;
            rhoe=@(s2) c2.*sqrt(1-s2.^2);
            S1e=@(s2,phi2) r1e(s2) +epps*rhoe(s2).*erho1(st,phi2);
            S2e=@(s2,phi2) r2e(s2) +epps*rhoe(s2).*erho2(st,phi2);
            S3e=@(s2,phi2) r3e(s2) +epps*rhoe(s2).*erho3(st,phi2);
            s3=@(s2) u2-st+s2;
            fxi2=fxi(st,pt);
            fyi2=fyi(st,pt);
            fzi2=fzi(st,pt);
            
            Cie=@(s2,phi2) 1./sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2);
            C11e=@(s2,phi2) (S1e(u2,pt)-S1e(s3(s2),phi2)).*(S1e(u2,pt)-S1e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            C12e=@(s2,phi2) (S1e(u2,pt)-S1e(s3(s2),phi2)).*(S2e(u2,pt)-S2e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            C13e=@(s2,phi2) (S1e(u2,pt)-S1e(s3(s2),phi2)).*(S3e(u2,pt)-S3e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            C22e=@(s2,phi2) (S2e(u2,pt)-S2e(s3(s2),phi2)).*(S2e(u2,pt)-S2e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            C23e=@(s2,phi2) (S2e(u2,pt)-S2e(s3(s2),phi2)).*(S3e(u2,pt)-S3e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            C33e=@(s2,phi2) (S3e(u2,pt)-S3e(s3(s2),phi2)).*(S3e(u2,pt)-S3e(s3(s2),phi2))./(sqrt((S1e(u2,pt)-S1e(s3(s2),phi2)).^2+(S2e(u2,pt)-S2e(s3(s2),phi2)).^2+(S3e(u2,pt)-S3e(s3(s2),phi2)).^2).^3);
            BI1e =@(x,y) (Cie(x,y)+C11e(x,y)).*fxi2 + C12e(x,y).*fyi2 + C13e(x,y).*fzi2;
            BI2e =@(x,y) (C12e(x,y)).*fxi2 + (Cie(x,y)+C22e(x,y)).*fyi2 + C23e(x,y).*fzi2;
            BI3e =@(x,y) (C13e(x,y)).*fxi2 + C23e(x,y).*fyi2 + (Cie(x,y)+C33e(x,y)).*fzi2;
            
            
%             if l==10
%                keyboard 
%             end
            
            if abs(st-u2)<0.0000001
                
                t11=quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-1,st,-pi,pt,'AbsTol',1e-3);
                t12=quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,1,-pi,pt,'AbsTol',1e-3);
                t13=quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-1,st,pt,pi,'AbsTol',1e-3);
                t14=quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,1,pt,pi,'AbsTol',1e-3);
                t21=quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-1,st,-pi,pt,'AbsTol',1e-3);
                t22=quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,1,-pi,pt,'AbsTol',1e-3);
                t23=quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-1,st,pt,pi,'AbsTol',1e-3);
                t24=quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,1,pt,pi,'AbsTol',1e-3);
                t31=quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-1,st,-pi,pt,'AbsTol',1e-3);
                t32=quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,1,-pi,pt,'AbsTol',1e-3);
                t33=quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-1,st,pt,pi,'AbsTol',1e-3);
                t34=quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,1,pt,pi,'AbsTol',1e-3);
                B1v(l+1)=t11+t12+t13+t14;
                B2v(l+1)=t21+t22+t23+t24;
                B3v(l+1)=t31+t32+t33+t34;
                
            elseif st-u2>0
                sig=+1;
                
                t11=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t12=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t13=sig*quad2d(@(x,y) BI1(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t14=sig*quad2d(@(x,y) -BI1e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t15=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t16=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t17=sig*quad2d(@(x,y) BI1(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t18=sig*quad2d(@(x,y) -BI1e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                t21=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t22=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t23=sig*quad2d(@(x,y) BI2(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t24=sig*quad2d(@(x,y) -BI2e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t25=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t26=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t27=sig*quad2d(@(x,y) BI2(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t28=sig*quad2d(@(x,y) -BI2e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                t31=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t32=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t33=sig*quad2d(@(x,y) BI3(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t34=sig*quad2d(@(x,y) -BI3e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t35=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t36=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t37=sig*quad2d(@(x,y) BI3(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t38=sig*quad2d(@(x,y) -BI3e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                
                B1v(l+1)=t11+t12+t13+t14+t15+t16+t17+t18;
                B2v(l+1)=t21+t22+t23+t24+t25+t26+t27+t28;
                B3v(l+1)=t31+t32+t33+t34+t35+t36+t37+t38;
            else
                sig=-1;
                
                t11=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t12=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t13=sig*quad2d(@(x,y) BI1(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t14=sig*quad2d(@(x,y) -BI1e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t15=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t16=sig*quad2d(@(x,y) BI1(x,y)-BI1e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t17=sig*quad2d(@(x,y) BI1(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t18=sig*quad2d(@(x,y) -BI1e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                t21=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t22=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t23=sig*quad2d(@(x,y) BI2(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t24=sig*quad2d(@(x,y) -BI2e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t25=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t26=sig*quad2d(@(x,y) BI2(x,y)-BI2e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t27=sig*quad2d(@(x,y) BI2(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t28=sig*quad2d(@(x,y) -BI2e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                t31=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-sig+st-u2,st,-pi,pt,'AbsTol',1e-3);
                t32=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,sig,-pi,pt,'AbsTol',1e-3);
                t33=sig*quad2d(@(x,y) BI3(x,y),-sig,-sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t34=sig*quad2d(@(x,y) -BI3e(x,y),sig,sig+st-u2,-pi,pt,'AbsTol',1e-3);
                t35=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),-sig+st-u2,st,pt,pi,'AbsTol',1e-3);
                t36=sig*quad2d(@(x,y) BI3(x,y)-BI3e(x,y),st,sig,pt,pi,'AbsTol',1e-3);
                t37=sig*quad2d(@(x,y) BI3(x,y),-sig,-sig+st-u2,pt,pi,'AbsTol',1e-3);
                t38=sig*quad2d(@(x,y) -BI3e(x,y),sig,sig+st-u2,pt,pi,'AbsTol',1e-3);
                
                
                B1v(l+1)=t11+t12+t13+t14+t15+t16+t17+t18;
                B2v(l+1)=t21+t22+t23+t24+t25+t26+t27+t28;
                B3v(l+1)=t31+t32+t33+t34+t35+t36+t37+t38;
                
            end
            
            %toc
            
        end
        
        s = linspace(-1+(ds/2),1-(ds/2),steps);
        phi=linspace(-pi+dphi/2,pi-dphi/2,steps+2);
        BIv=zeros(3,steps+2,steps);
        for l=0:((steps+2)*steps-1)
            j=floor(l/(steps+2))+1;
            i=rem(l,steps+2)+1;
            fv=tv(i,j,1)*fxi(S(i,j),Phi(i,j)) +tv(i,j,2)*fyi(S(i,j),Phi(i,j))+tv(i,j,3)*fzi(S(i,j),Phi(i,j));
            BIv(1,i,j)=B1v(l+1)+(Z2(a(i,j),alpha(i,j))*fxi(S(i,j),Phi(i,j)) + tv(i,j,1)*fv*(Z1(a(i,j),alpha(i,j))-Z2(a(i,j),alpha(i,j))));
            BIv(2,i,j)=B2v(l+1)+(Z2(a(i,j),alpha(i,j))*fyi(S(i,j),Phi(i,j)) + tv(i,j,2)*fv*(Z1(a(i,j),alpha(i,j))-Z2(a(i,j),alpha(i,j))));
            BIv(3,i,j)=B3v(l+1)+(Z2(a(i,j),alpha(i,j))*fzi(S(i,j),Phi(i,j)) + tv(i,j,3)*fv*(Z1(a(i,j),alpha(i,j))-Z2(a(i,j),alpha(i,j))));
        end
        
        q1=reshape(BIv(1,:,:),steps+2,steps)-q1m;
        q2=reshape(BIv(2,:,:),steps+2,steps)-q2m;
        q3=reshape(BIv(3,:,:),steps+2,steps)-q3m;
        
        Mipal=zeros(3,steps+2,steps);
        tq=tv(:,:,1).*q1+tv(:,:,2).*q2+tv(:,:,3).*q3;
        Mipal(1,:,:)=(q1./z2p +tv(:,:,1).*tq.*(1./z1p -1./z2p));
        Mipal(2,:,:)=(q2./z2p +tv(:,:,2).*tq.*(1./z1p -1./z2p));
        Mipal(3,:,:)=(q3./z2p +tv(:,:,3).*tq.*(1./z1p -1./z2p));
        
        Q= reshape(dphi*sum(Mipal,2),3*steps,1);
        
        A0=Bl*Q;
        
        fb=B\A0;
        
        MT=reshape(A0- Bl*fb,3,steps);
        
        fi = zeros(3,steps+2,steps);
        li = zeros(3,steps+2,steps);
        for i=1:(steps+2)
            for j=1:(steps)
                fi(:,i,j)= Mipal(:,i,j) -(t(s(j)).*t(s(j))')*MT(:,j)/z1p(i,j) -((eye(3)-t(s(j)).*t(s(j))'))*MT(:,j)/z2p(i,j);
                li(:,i,j)= cross([S1(S(i,j),Phi(i,j)),S2(S(i,j),Phi(i,j)),S3(S(i,j),Phi(i,j))],fi(:,i,j));
            end
        end
        
        f = f+((-1)^n)*fi;
        F(1,n+1)=F(1,n) + ((-1)^n)*sum(fi(1,:,:),'all')*ds*dphi;
        F(2,n+1)=F(2,n) + ((-1)^n)*sum(fi(2,:,:),'all')*ds*dphi;
        F(3,n+1)=F(3,n) + ((-1)^n)*sum(fi(3,:,:),'all')*ds*dphi;
        L(1,n+1)=L(1,n) + ((-1)^n)*sum(li(1,:,:),'all')*ds*dphi;
        L(2,n+1)=L(2,n) + ((-1)^n)*sum(li(2,:,:),'all')*ds*dphi;
        L(3,n+1)=L(3,n) + ((-1)^n)*sum(li(3,:,:),'all')*ds*dphi;
        
        fxi =@(y,x) interp2(phi,s,reshape(fi(1,:,:),steps+2,steps)',x,y,'spline');
        fyi =@(y,x) interp2(phi,s,reshape(fi(2,:,:),steps+2,steps)',x,y,'spline');
        fzi =@(y,x) interp2(phi,s,reshape(fi(3,:,:),steps+2,steps)',x,y,'spline');
        q1m=q1;
        q2m=q2;
        q3m=q3;
    end
end


end

