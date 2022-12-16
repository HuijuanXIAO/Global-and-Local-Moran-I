%Matlab code for Global Bivariate Moran Index
clear all;
A=csvread('Data.csv',1,0);
y=A(:,[1]); 
x=A(:,[2]); 
W=csvread('Weight.csv',1,0);
w=W;
n=length(y);
xs=std(x,1); 
ys=std(y,1); 
xm=mean(x);
ym=mean(y);

%normalizes each row.
for k =1:n
    w(k,:) = w(k,:)/norm(w(k,:),1);
end

%Global Moran's I
g=0;
a=0;
for i=1:1:n
    for j=1:1:n
         if (i==j)
            continue  
        end
        g=g+w(i,j)*(x(i)-xm)*(y(j)-ym);
        a=a+w(i,j);                   
    end
end
Global=g/(xs*ys*a);
k21=0;
k22=0;
s0=0;
s1=0;
s2=0;
for i=1:1:n
   k21=k21+(x(i)-xm)^2*(y(i)-ym)^2;
   k22=k22+(x(i)-xm)*(y(i)-ym);
     for j=1:1:n
        s0=s0+w(i,j);
        s1=s1+(w(i,j)+w(j,i))^2;
     end
end
for i=1:1:n
    s12(i)=0;
    s21(i)=0;
  for j=1:1:n
     s12(i)=s12(i)+w(i,j);
     s21(i)=s21(i)+w(j,i);
  end
    s2=s2+(s12(i)+s21(i))^2;
end
b2=(n*k21)/(k22^2);

s1=s1/2;
ei=-1/(n-1);  
vari_G=(n*((n^2-3*n+3)*s1-n*s2+3*s0^2)-b2*((n^2-n)*s1-2*n*s2+6*s0^2))/(s0^2*(n-1)*(n-2)*(n-3))-ei^2; 
z_test_G=(Global-ei)/vari_G^(1/2); 
p_value_G = 2 *(1- normcdf(abs(z_test_G)));