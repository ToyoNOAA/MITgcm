% This is a matlab script that generates the input data

% Dimensions of grid
nx=200.0; ny=200.0; nz=50.0;
% Nominal depth of model (meters)
H=1000.0;
% Size of domain (m)
Lx=8.0e4; Ly=8.0e4;
% Resolution (m)
dx=Lx/nx; dy=Ly/ny; dz=H/nz;
% Rotation
f=1.e-4;
% Stratification
N=0.0;
% Surface temperature
Ts=20.;
% Lapse rate
dTdz=-0.01;
% Gravity
g=10.;
% E.O.S.
alpha=2.e-4; Tz=N^2/(g*alpha);

% Set vortex radius, vertical decay scale and max speed
R=1.0e4; Ht=250; Um=0.1;

% Compute vortex horizontal decay scale, and velocity and temperature amplitudes
Lt=2*R^2;
A=sqrt(exp(1))*Um/R;
t=A*f*Lt/(2*g*alpha*Ht);

% Print vertical levels for /input/data
sprintf('delZ = %d * %7.6g,',nz,dz)

% Initial temperature and velocity distributions
x=zeros(nx,1); y=zeros(ny,1); z=zeros(nz,1);

for i=1:nx
        x(i)=(i-1)*dx;
end
for i=1:ny
	y(i)=(i-1)*dy;
end
for i=1:nz
	z(i)=(i-1)*dz;
end

T=zeros(nx,ny,nz); U=zeros(nx,ny,nz); V=zeros(nx,ny,nz);

for k=1:nz
	for j=1:ny
		for i=1:nx
			r(i,j) = norm([x(i)-Lx/2 y(j)-Ly/2]);
			T(i,j,k)=Ts+dTdz*z(k); %background temperature
			T(i,j,k)=T(i,j,k)+t*exp(-z(k)/Ht)*exp(-r(i,j)^2/Lt); %perturbed temperature
			U(i,j,k)=-A*y(j)*exp(-z(k)/Ht)*exp(-r(i,j)^2/Lt);
			V(i,j,k)=A*x(i)*exp(-z(k)/Ht)*exp(-r(i,j)^2/Lt);
		end
	end
end

fid=fopen('T.bin','w','b'); fwrite(fid,T,'real*4'); fclose(fid);
fid=fopen('U.bin','w','b'); fwrite(fid,U,'real*4'); fclose(fid);
fid=fopen('V.bin','w','b'); fwrite(fid,V,'real*4'); fclose(fid);
