function [ao,a,b,wn]= fourier_series_num(X,Y)
f_bins=500;                         % no. of freq. bins

%X=0:1e-3:1;                   % time vector
%Y=sin(2*pi*10*X);        % change the function as u like (or put data points)
n=length(X);                 % no. of samples
f_s=1/(X(5)-X(4));                 % sampling frequency

ao=(1/(X(end)-X(1)))*trapz(Y)*(1/f_s);         %fourier coeff. ao , by the integration int(f(x),0,L)*(1/L)

cc=zeros(f_bins,n);
ss=zeros(f_bins,n);

for nn=1:f_bins
   for z=1:length(X);
        zz=X(z);
         cc(nn,z)=cos(pi*nn*zz/(X(end)-X(1)));
         ss(nn,z)=sin(pi*nn*zz/(X(end)-X(1))); 
    end
end

a=zeros(f_bins,1);                   
b=zeros(f_bins,1);                  
aa=zeros(f_bins,n);
bb=zeros(f_bins,n);

for ii=1:f_bins
    a(ii)=(1/(X(end)-X(1)))*trapz(Y.*cc(ii,:))*(1/f_s);  % fourier coeff. an
    b(ii)=(1/(X(end)-X(1)))*trapz(Y.*ss(ii,:))*(1/f_s);  % fourier coeff. bn
     for z=1:length(X);
       zz=X(z);
       aa(ii,z)=a(ii)*cos(pi*ii*zz/(X(end)-X(1)));
       bb(ii,z)=b(ii)*sin(pi*ii*zz/(X(end)-X(1)));
     end
end

f=.5*ao+sum(aa,1)+sum(bb,1);                    % fourier series
amp=sqrt(a.^2+b.^2);
frequency=(1:f_bins)/(X(end)-X(1));                      % corresponding frequency values 
wn=2*pi*frequency;
amp=[ao;amp(1:end)]';                        % amplitude spectrum

figure
plot(X,Y)
hold on
plot(X,f,'s')
legend('Real signal','Fourier representation')
