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

% Set vortex radius, vertical decay scale and maximum velocity
R=Lx/8; Ht=250; Um=0.1;

% Compute correction term, vorticity and temperature amplitude
C=R^2*(log(R^2)-1); w=-2*Um/R; t=-w*f/(4*g*alpha*Ht);

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
			r(i,j)=norm([x(i)-Lx/2 y(j)-Ly/2]);
			T(i,j,k)=Ts+dTdz*z(k);
			if r(i,j) < R
				U(i,j,k)=-(w/2)*(y(j)-Ly/2)*exp(-z(k)/Ht);
				V(i,j,k)=(w/2)*(x(i)-Lx/2)*exp(-z(k)/Ht);
				T(i,j,k)=T(i,j,k)+t*(r(i,j)^2+C)*exp(-z(k)/Ht);
			else
				U(i,j,k)=-(w*R^2/2)*(y(j)-Ly/2)*exp(-z(k)/Ht)/r(i,j)^2;
				V(i,j,k)=(w*R^2/2)*(x(i)-Lx/2)*exp(-z(k)/Ht)/r(i,j)^2;
				T(i,j,k)=T(i,j,k)+t*R^2*log(r(i,j)^2)*exp(-z(k)/Ht);
			end
		end
	end
end

fid=fopen('T.bin','w','b'); fwrite(fid,T,'real*4'); fclose(fid);
fid=fopen('U.bin','w','b'); fwrite(fid,U,'real*4'); fclose(fid);
fid=fopen('V.bin','w','b'); fwrite(fid,V,'real*4'); fclose(fid);
