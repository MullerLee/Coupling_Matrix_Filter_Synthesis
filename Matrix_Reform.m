% filename: Matrix_Reform.m
% author:   Zve.L
% date:     3/13/2020 
% rev.log   
% Please Run After Coupling_Matrix_Extraction.m

%Process of Coupling Matrix M
%        R2-------R4
%            R3
%   S----R1-------R5-----L
%   0       M12     M13     0       M15 
%   M21     0       M23     M24     0
%   M31     M32     0       M34     M35
%   0       M42     M43     0       M45
%   M51     0       M53     M54     0

M0=[0.1,0.4,0.2,0.3,0.5;
   0.4,0.7,0.2,0.5,0.8;
   0.2,0.2,0.8,0.3,0.4;
   0.3,0.5,0.3,0.7,0.2;
   0.5,0.8,0.4,0.2,0.1;];
M=M0;
% Operate M14, pivot [3,4] 
k2=1; l2=4; m2=1; n2=3; c2=-1; %(in row)
theta = atan(c2*M(k2,l2)/M(m2,n2));
cr=cos(theta);
sr=sin(theta);
MM=M;
MM(3,1)=cr*M(3,1)-sr*M(4,1); MM(3,2)=cr*M(3,2)-sr*M(4,2); MM(3,5)=cr*M(3,5)-sr*M(4,5);
MM(4,1)=sr*M(3,1)+cr*M(4,1); MM(4,2)=sr*M(3,2)+cr*M(4,2); MM(4,5)=sr*M(3,5)+cr*M(4,5);
MM(1,3)=cr*M(1,3)-sr*M(1,4); MM(2,3)=cr*M(2,3)-sr*M(2,4); MM(5,3)=cr*M(5,3)-sr*M(5,4);
MM(1,4)=sr*M(1,3)+cr*M(1,4); MM(2,4)=sr*M(2,3)+cr*M(2,4); MM(5,4)=sr*M(5,3)+cr*M(5,4);
M=MM;

% Operate M25, pivot [3,5] 
k2=2; l2=5; m2=2; n2=3; c=-1;%(in row)
theta = atan(c2*M(k2,l2)/M(m2,n2));
cr=cos(theta);
sr=sin(theta);
MM=M;
MM(3,1)=cr*M(3,1)-sr*M(5,1); MM(3,2)=cr*M(3,2)-sr*M(5,2); MM(3,4)=cr*M(3,4)-sr*M(5,4);
MM(5,1)=sr*M(3,1)+cr*M(5,1); MM(5,2)=sr*M(3,2)+cr*M(5,2); MM(5,4)=sr*M(3,4)+cr*M(5,4);
MM(1,3)=cr*M(1,3)-sr*M(1,5); MM(2,3)=cr*M(2,3)-sr*M(2,5); MM(4,3)=cr*M(4,3)-sr*M(4,5);
MM(1,5)=sr*M(1,3)+cr*M(1,5); MM(2,5)=sr*M(2,3)+cr*M(2,5); MM(4,5)=sr*M(4,3)+cr*M(4,5);
M=MM;