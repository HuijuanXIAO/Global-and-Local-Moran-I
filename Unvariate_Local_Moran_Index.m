%Matlab code for Local Unvariate Moran Index
clear all;
A=csvread('Data.csv',1,0);
x=A(:,[1]); 
W=csvread('Weight.csv',1,0);
w=W;
n=length(x);
s=var(x,1); 
m=mean(x);
y=0;
a=0;
%normalizes each row.
for k =1:n
    w(k,:) = w(k,:)/norm(w(k,:),1);
end

a=0;
for i=1:1:n
    for j=1:1:n
         if (i==j)
            continue  
         end
         a=a+w(i,j);                
    end
end

%Local Moran's I
L=zeros([length(y),1]);
for i=1:n
    L_tmp=0;
    for j=1:n
        L_tmp=L_tmp+w(i,j)*(x(i)-m)*(x(j)-m);
    end
    L(i,1)=L_tmp;
end

Local=n*L/(s*a); 

k21=0;
k22=0;
s0=0;
s1=0;
s2=0;
for i=1:1:n
   k21=k21+(x(i)-m)^4;
   k22=k22+(x(i)-m)^2;
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

Wi_2=zeros([n,1]);
for i=1:n
    Wi_2_tmp=0;
    for j=1:n
        if (j==i)
            continue 
        end
        Wi_2_tmp=Wi_2_tmp+w(i,j)^2;
    end
    Wi_2(i,1)=Wi_2_tmp;
end

Wi=zeros([n,1]);
for i=1:n
    Wi_tmp=0;
    for j=1:n
        if (j==i)
            continue  
        end
        Wi_tmp=Wi_tmp+w(i,j);
    end
    Wi(i,1)=Wi_tmp;
end

ei_L=zeros(n,1); 
for i=1:n
    ei_L(i,1)=-Wi(i,1)./(n-1);  
end

vari_L=zeros(n,1);
for i=1:n
    vari_L(i,1)=Wi_2(i,1).*(n-b2)./(n-1)+(((Wi(i,1))^2-Wi_2(i,1))*(2*b2-n))./((n-1)*(n-2))-(ei_L(i,1)).^2;
end

z_test_L=zeros([n,1]); 
for i=1:n
    z_test_L(i,1)=(Local(i,1)-ei_L(i,1))/(vari_L(i,1).^(1/2)); 
end

p_value_L=zeros([n,1]); 
for i=1:n
    p_value_L(i,1) = 2 *(1- normcdf(abs(z_test_L(i,1))));
end
