function [q,L,T,angles,res,so,Cqq,Cpar,xyz_new]=DQA_3D_2D_transformation(xyz1,xyz2,P)
%
% DQA_3d_transformation (Asymmetric similarity) with dual quaternions
% using scaled quaternions
% reverse and direct problem

% Model: xyz2=(2*w_r'*s)+(w_r'*Q_r*xyz1)  
% function performs both reverse and direct problem-solving based on dual quaternions between two 3D Cartesian coordinate system
% The function can also 2D transformation between given coordinate systems.
% INPUT:
% n :number of control points
% xyz1 = nx3 start system coordinate matrix ([x1 y1 z1])[in meter]
% xyz1 contain [ n control points xyz1]
%       | x1 y1 z1 |
%  xyz1=| x2 y2 z2 |
%       | ........ |
%       | xn yn zn |
% xyz2 = nx3 target system coordinate matrix ([x2 y2 z2])[in meter]
% xyz2 contain [ n control points xyz2 ]
%       | X1 Y1 Z1 |
%  xyz2=| X2 Y2 Z2 |
%       | ........ |
%       | Xn Yn Zn |
%  if there are new point to be transformed (revers problem+direct problem)
% xyz1 = mx3 start system coordinate matrix
% the new points are added under the control points in the xyz1 matrix.
% m :number of total points (control points+new points to be transformed)
%       | x1    y1    z1 |==>control points
%       | x2    y2    z2 |==>control points
%       | .............. |==>control points
%  xyz1=| xn    yn    zn |==>control points
%       | xn+1 yn+1 zn+1 |==>new points  
%       | .............. |==>new points
%       | xm   ym   zm   |==>new points
% (m-n)= new points to be transformed (if any) ]
% FOR 2D TRANSFORMATIONS
% All of Z1 and Z2 coordinates must be entered as zero.
% P = nx1 weighted of points 
% m-n : number of new points to be transformed
% m must be greater than or equal n (m>=n)
% if m=n there is no new point to be transformed (revers problem)
% if m>n there is new point to be transformed (revers problem+direct problem)
%  OUTPUT:
% q =[2 x 8] matrix  dual quaternion elements [r0 r1 r2 r3 s0 s1 s2 s3] [unitless]
%       first row contained scaled quaternions
%       second row contained unit quaternions
% L = scale factor [unitless]
% T = translation vector (Tx,Ty,Tz)[in meter]
% angles = (e,p,w) rotation angles (x,y,z) [in radian]
% res = residual of coordinates [in meter]
% so = standard error of transformation [in meter]
% Cqq = covariance matrix of scaled quaternions
% (r0,r1,r2,r3,s0,s1,s2,s3)respectively
% Cpar = covariance matrix of transformation parameters
% (L,e,p,w,tx,ty,tz,0)respectively
%
% xyz_new=(m-n)x3 matrix is the coordinates of the new points calculated in the target system.  (if m>n)
% xyz_new=[]  (if m=n)
%==========================================================================
%==========================================================================
% SAMPLE INPUT DATA For 3D Transformations
% n=4 number of control point
% m=5 number of row of xyz1 
% m-n=1 number of new point to be converted (if any)
% xyz1=[ 6432.5800 7254.1200 200.6000 
%        6354.3700 5724.5800 174.5700 
%        7221.4400 6355.0800 254.5800 
%        6433.5800 7255.1200 201.6000 
%        6533.5800 7055.1200 211.6000 ] this row belong to new points to be converted.
% xyz2=[ 4208.8321 2111.9343 4182.9434 
%        2034.5929 2073.9091 4924.8221 
%        3397.0341 1919.6811 5773.1190 
%        4207.7300 2110.9400 4182.1100 ] 
%
% P=[1 2 3 4] weight
%======================================
% SAMPLE INPUT DATA For 2D Transformations
% n=4 number of control point
% m=5 number of row of xyz1 
% m-n=1 number of new point to be converted (if any)
%  xyz1=[  9043.7400    5208.7900       0      
%    	   9218.4200    4833.4900       0   
%    	   9000.0000    5000.0000       0   
%    	   9220.0200    5166.9100       0   
%          9242.7000    5039.3800       0   
%          9106.1700    5050.7100       0] this row belong to new points to be % converted.
%  xyz2= [ 4618.7200    4068.8300       0
%          5579.4100    1115.6000       0
%          4103.9800    2553.3800       0   
%    	   5893.3800    3597.0300       0   
%    	   5946.7000    2626.7000       0 ]
% P=[1 1 1 1] weight
% ======================================================
% Please refer to:
% 
% Bektas, S.(2023) An expanded dual quaternion algorithm for 3D Helmert %transformation and determination of the VCV matrix of the transformation’s %parameters. Journal of Spatial Science, DOI: 10.1080/14498596.2023.2274997
% Bektas, S.(2017) Adjustment Calculus,Ondokuz Mayis University press,
% ISBN 975-7636-54-1 Samsun,Turkey 
%
%==========================================================================
% written by Sebahattin BEKTAS January 2023
% sbektas@omu.edu.tr
   tol=10^-8;
   ro = 180 / pi;
     [n,sut]=size(xyz2);
     [m,sut]=size(xyz1);
   
ons=n;% control(common) point number to be converted
yns = m-n; % new point number to be converted
      
    X1=xyz1(:,1);
    Y1=xyz1(:,2);
    Z1=xyz1(:,3);
    X2=xyz2(:,1);
    Y2=xyz2(:,2);
    Z2=xyz2(:,3);
    pp=P;
      
P=eye(3*ons+1);
for k=1:ons
    kk=3*k;
    
    P(kk-2,kk-2)=pp(k);
    P(kk-1,kk-1)=pp(k);
    P(kk-0,kk-0)=pp(k);
end
%P
dimen=3; %dimension of transformation
np=7; % number of transfomation parameters
if(sum(abs(Z2))<0.01)
dimen=2;np=4;
end
%fprintf(dimen,'D TRANSFORMATION\n', %5d
fprintf('%10d D COORDINATE TRANSFORMATION!!! ATTENTION!!! \n',dimen)
%
  xx=xyz1(:,1);yy=xyz1(:,2);zz=xyz1(:,3);%first system
  x_new=xyz2(:,1);y_new=xyz2(:,2);z_new=xyz2(:,3); %second system
  
  %  fprintf ('initial values of  quaternions\n')
   
  r0=1  ;r1=-0.;r2=-0.;r3=-0.;   d0=0  ;d1=-0.;d2=-0.;d3=-0.;
  
vek=[r0,r1,r2,r3,d0,d1,d2,d3];
it=0;did=0;L=0;
while did==0
    it=it+1;
       for i=1:ons
       ib=3*(i-1)+1;
       x=xx(i);y=yy(i);z=zz(i);
                                                                                                                                                                                                                                    
 pu =[ (z)*(2*r0*r2 + 2*r1*r3) - (y)*(2*r0*r3 - 2*r1*r2) + 2*d0*(r0) - 2*d1*(r3) + 2*d2*(r2) - 2*d3*(r1) + (x)*(r0^2 + r1^2 - r2^2 - r3^2)
 (x)*(2*r0*r3 + 2*r1*r2) - (z)*(2*r0*r1 - 2*r2*r3) + 2*d1*(r0) + 2*d0*(r3) - 2*d2*(r1) - 2*d3*(r2) + (y)*(r0^2 - r1^2 + r2^2 - r3^2)
 (y)*(2*r0*r1 + 2*r2*r3) - (x)*(2*r0*r2 - 2*r1*r3) - 2*d0*(r2) + 2*d1*(r1) + 2*d2*(r0) - 2*d3*(r3) + (z)*(r0^2 - r1^2 - r2^2 + r3^2)
                                                                                                   2*d0*(r1) + 2*d1*(r2) + 2*d3*(r0) + 2*d2*(r3)];
 
 
Qjacob =[ 2*d0 + 2*r0*(x) - 2*r3*(y) + 2*r2*(z), 2*r1*(x) - 2*d3 + 2*r2*(y) + 2*r3*(z), 2*d2 - 2*r2*(x) + 2*r1*(y) + 2*r0*(z), 2*r1*(z) - 2*r3*(x) - 2*r0*(y) - 2*d1,  2*(r0), -2*(r3),  2*(r2), -2*(r1)
 2*d1 + 2*r3*(x) + 2*r0*(y) - 2*r1*(z), 2*r2*(x) - 2*d2 - 2*r1*(y) - 2*r0*(z), 2*r1*(x) - 2*d3 + 2*r2*(y) + 2*r3*(z), 2*d0 + 2*r0*(x) - 2*r3*(y) + 2*r2*(z),  2*(r3),  2*(r0), -2*(r1), -2*(r2)
 2*d2 - 2*r2*(x) + 2*r1*(y) + 2*r0*(z), 2*d1 + 2*r3*(x) + 2*r0*(y) - 2*r1*(z), 2*r3*(y) - 2*r0*(x) - 2*d0 - 2*r2*(z), 2*r1*(x) - 2*d3 + 2*r2*(y) + 2*r3*(z), -2*(r2),  2*(r1),  2*(r0), -2*(r3)];
%                                              2*d3,                                              2*d0,                                              2*d1,                                              2*d2,  2*(r1),  2*(r2),  2*(r3),  2*(r0)];
                 
A(ib:ib+2,:)=Qjacob;
sab(ib+0)=x_new(i)-pu(1);
sab(ib+1)=y_new(i)-pu(2);
sab(ib+2)=z_new(i)-pu(3);
%sab(ib+3)=0-pu(4);
 
 
       end
   A(3*ons+1,:)    =[ 2*d3    2*d0  2*d1  2*d2  2*(r1)  2*(r2)  2*(r3)  2*(r0)];
   sab(3*ons+1)= -(2*d0*(r1) + 2*d1*(r2) + 2*d3*(r0) + 2*d2*(r3));
             %  A(4*ons+1,:)=[d3  d0  d1  d2   r1  r2  r3  r0 ]; sab(4*ons+1)=(-d0*r1-d1*r2-d2*r3-d3*r0);  
    % A,sab
       L=(r0^2 + r1^2 + r2^2 + r3^2);
  
   d_vek=(sqrt(P)*A)\(sqrt(P)*sab');
r0=r0+d_vek(1);
r1=r1+d_vek(2);
r2=r2+d_vek(3);
r3=r3+d_vek(4);
d0=d0+d_vek(5);
d1=d1+d_vek(6);
d2=d2+d_vek(7);
d3=d3+d_vek(8);
L=(r0^2 + r1^2 + r2^2 + r3^2);
vek=[r0,r1,r2,r3,d0,d1,d2,d3,L];
 
normx=norm(d_vek);
if it>=30| normx<tol
    did=1;
end
end
 
   w_r =[ r0        r3       -r2 r1
         -r3        r0        r1 r2
          r2       -r1        r0 r3
         -r1       -r2       -r3 r0];
  d=[d0,d1,d2,d3]'; T_ilk=2*w_r'*d;
    q_scaled=[r0,r1,r2,r3,d0,d1,d2,d3];
       r0=r0/sqrt(L);r1=r1/sqrt(L);r2=r2/sqrt(L);r3=r3/sqrt(L);   d0=d0*sqrt(L);d1=d1*sqrt(L);d2=d2*sqrt(L);d3=d3*sqrt(L);
       q_unit=[r0,r1,r2,r3,d0,d1,d2,d3];
       q=[q_scaled;q_unit];
         pu =[ (z)*(2*r0*r2 + 2*r1*r3) - (y)*(2*r0*r3 - 2*r1*r2) + 2*d0*(r0) - 2*d1*(r3) + 2*d2*(r2) - 2*d3*(r1) + (x)*(r0^2 + r1^2 - r2^2 - r3^2)
 (x)*(2*r0*r3 + 2*r1*r2) - (z)*(2*r0*r1 - 2*r2*r3) + 2*d1*(r0) + 2*d0*(r3) - 2*d2*(r1) - 2*d3*(r2) + (y)*(r0^2 - r1^2 + r2^2 - r3^2)
 (y)*(2*r0*r1 + 2*r2*r3) - (x)*(2*r0*r2 - 2*r1*r3) - 2*d0*(r2) + 2*d1*(r1) + 2*d2*(r0) - 2*d3*(r3) + (z)*(r0^2 - r1^2 - r2^2 + r3^2)
                                                                                                   2*d0*(r1) + 2*d1*(r2) + 2*d3*(r0) + 2*d2*(r3)];
 
  
      
d=[d0,d1,d2,d3]';
s=[r0,r1,r2,r3]';
w_r =[ r0        r3       -r2 r1
         -r3        r0        r1 r2
          r2       -r1        r0 r3
         -r1       -r2       -r3 r0];
     
 wtq=[       r0*(r0) + r1*(r1) - r2*(r2) - r3*(r3),                   2*r1*(r2) - r0*(r3) - r3*(r0),                   r0*(r2) + r2*(r0) + 2*r1*(r3),       r1*(r0) - r0*r1 - r2*(r3) + r3*(r2)
                   r0*(r3) + 2*r2*(r1) + r3*(r0),            r0*(r0) - r1*(r1) + r2*(r2) - r3*(r3),                   2*r2*(r3) - r1*(r0) - r0*(r1),       r2*(r0) - r0*r2 + r1*(r3) - r3*(r1)
                   2*r3*(r1) - r2*(r0) - r0*(r2),                   r0*(r1) + r1*(r0) + 2*r3*(r2),       r0*(r0) - r1*(r1) - r2*(r2) + r3*(r3),       r2*(r1) - r1*(r2) - r0*r3 + r3*(r0)
 r0*(r1) - (r0)*(r1) - r2*(r3) + r3*(r2), r0*(r2) - (r0)*(r2) + r1*(r3) - r3*(r1), r0*(r3) - (r0)*(r3) - r1*(r2) + r2*(r1), r0*(r0) + r1*(r1) + r2*(r2) + r3*(r3)];
    
fprintf('Transformed Control Points Coordinates\n')
   for i=1:ons
       ib=4*(i-1)+1;
       x=xx(i);y=yy(i);z=zz(i);
       pu=   (2*w_r'*d)+L*(wtq*[x;y;z;0]);
 
  fprintf('%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f  %8.1f\n',xx(i),yy(i),zz(i), x_new(i),y_new(i),z_new(i),pu)
   end
     if m>n
          fprintf ('Transformed New Points Coordinates \n')
   for j=ons+1:m
       i=j-ons;
       ib=4*(i-1)+1;
    %   x=xx(i);y=yy(i);z=zz(i);
       x=xx(j);y=yy(j);z=zz(j);
   
    
  pu=   (2*w_r'*d)+L*(wtq*[x;y;z;0]);
       
    fprintf('%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %8.1f\n',xx(j),yy(j),zz(j), pu)
    xyz_new(i,:)=pu;
   end
end
if m==n
    xyz_new=[];
end
T=2*w_r'*d;
ax = atan2 (-2*(r0*r1+r2*r3),(r0^2-r1^2-r2^2+r3^2));
ay = asin (2*(-r0*r2+r3*r1));
az = atan2 (-2*(r0*r3+r2*r1),(r0^2+r1^2-r2^2-r3^2)) ;
angles=[ax,ay,az];
    
    V=A*d_vek+sab';
    
    ij=1;
    for i=1:ons
        res(i,1)=V(ij);ij=ij+1;
        res(i,2)=V(ij);ij=ij+1;
        res(i,3)=V(ij);ij=ij+1;
    end
    
 vtpv=V'*P*V;
so=sqrt(vtpv/(dimen*ons-np));
Q=(A'*(P)*A)^-1;
Cqq=so^2*Q;
aci=[ax,ay,az]*180/pi;
cx=cos(ax);sx=sin(ax);
cy=cos(ay);sy=sin(ay);
cz=cos(az);sz=sin(az);
q0=r0*sqrt(L); q1=r1*sqrt(L); q2=r2*sqrt(L); q3=r3*sqrt(L);
d0=d0/sqrt(L); d1=d1/sqrt(L); d2=d2/sqrt(L); d3=d3/sqrt(L);
Fq=2*[q0 q1 q2 q3
    -(q1*cx+q0*sx)/(L*cy) -(q0*cx-q1*sx)/(L*cy) -(q3*cx-q2*sx)/(L*cy) -(q2*cx+q3*sx)/(L*cy)
    -(q2+q0*sy)/(L*cy) -(-1*q3+q1*sy)/(L*cy) -(q0+q2*sy)/(L*cy) -(-1*q1+q3*sy)/(L*cy)
    -(q3*cz+q0*sz)/(L*cy) -(q2*cz+q1*sz)/(L*cy) -(q1*cz-q2*sz)/(L*cy) -(q0*cz-q3*sz)/(L*cy)];
Fd =[ 2*d0, -2*d3,  2*d2, -2*d1,  2*(q0), -2*(q3),  2*(q2), -2*(q1)
 2*d1, -2*d2, -2*d3,  2*d0,  2*(q3),  2*(q0), -2*(q1), -2*(q2)
 2*d2,  2*d1, -2*d0, -2*d3, -2*(q2),  2*(q1),  2*(q0), -2*(q3)
 2*d3,  2*d0,  2*d1,  2*d2,  2*(q1),  2*(q2),  2*(q3),  2*(q0)];
%%%%%%%%%%%%% covariance matrix
Ft=[Fq  zeros(4,4);  Fd];
Cpar=Ft*Cqq*Ft';
  end