%Appendix: Implementation of asynchronous 2D-TDOA localization algorithm
%for simulations in MATLAB.
%In this appendix, the tested implementation in Matlab of our 2D-TDOA 
%localization algorithm is given for the easier repetition of the obtained 
%results and the future hardware implementation, due to the complexity of the formulas (25)-(31) in the paper.

function [res] = lokacijaTDOA (a,b,c,d,e) 
% Outpout: res
% where: res(1)=x* & res(2) = y*;
% x*, y* solution of 2D-TDOA shown in Fig. 3 of the paper.
% Inputs: a, b, c, d, e -> for its meaning look in the paper. 
% Defaults for x = 3.0, y=4.0 are (described in the paper):
% a = 2.211102550927978 ; 
% b = 7.369316876852981; 
% c = 9; d = 15; e = 7;
A=4*(e^2-b^2)*(c^2-a^2)/(a^2)+4*(d^2-b^2); 
B=4*(a*a*c-c*c*c)*(e^2-b^2)/(a*a)+4*(d*b*b-d^3-d*e*e); 
C=b^4-2*b*b*d*d+2*d*d*e*e+d^4+e^4+(e^2-b^2)*(a^4+c^4-2*a*a*c*c)/(a*a)-2*b*b*e*e; 
D=64*d*d*e*e*(c*c-a*a)/(a*a); 
E=(64*d*d*e*e*(a*a*c-c*c*c)+64*d*e*(c*c-a*a)*(b*b*e-e*d*d-e^3))/(a*a); 
F=(16*((e^3+e*d*d-e*b*b)^2)*(c^2-a^2)+64*d*e*(b*b*e-e*d*d-e^3)*(a*a*c-c^3)+16*d*d*e*e*(a^4+c^4-2*a*a*c*c))/(a*a); 
G=(16*((e^3+e*d*d-e*b*b)^2)*(a*a*c-c^3)+16*d*e*(b*b*e-e*d*d-e^3)*(a^4+c^4-2*a*a*c*c))/(a*a); 
H=4*((a^4+c^4-2*a*a*c*c)*((e^3+e*d*d-e*b*b)^2))/(a*a); 
AA=A^2-D; 
BB=2*A*B-E; 
CC=2*A*C+B^2-F; 
DD=2*B*C-G; 
EE=C^2-H;  
u=EE/AA;
p=BB/AA;
q=CC/AA; 
r=DD/AA;  
x1 = -(1/2) *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3) - (1/2) *sqrt(-(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) - (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/2 - (-p^3 + 4 *p *q - 8 *r)/(4 *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3)) - (4 *q)/3) - p/4;
x2 = -(1/2)* sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3) + (1/2) *sqrt(-(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) - (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/2 - (-p^3 + 4 *p *q - 8 *r)/(4 *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3)) - (4 *q)/3) - p/4;
x3 = (1/2) *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3) - (1/2) *sqrt(-(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) - (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/2 + (-p^3 + 4 *p *q - 8 *r)/(4 *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3)) - (4 *q)/3) - p/4; 
x4 = (1/2) *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3) + (1/2) *sqrt(-(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) - (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/2 + (-p^3 + 4 *p *q - 8 *r)/(4 *sqrt((sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)/(3 *2^(1/3)) + (2^(1/3) *(-3 *p *r + q^2 + 12 *u))/(3 *(sqrt((27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^2 - 4 *(-3 *p *r + q^2 + 12 *u)^3) + 27 *p^2 *u - 9 *p *q *r + 2 *q^3 - 72 *q *u + 27 *r^2)^(1/3)) + p^2/4 - (2 *q)/3)) - (4 *q)/3) - p/4; 
y1=(1/(2*a))*sqrt(4*x1^2*(c^2-a^2)+4*x1*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2);
y2=(1/(2*a))*sqrt(4*x2^2*(c^2-a^2)+4*x2*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
y3=(1/(2*a))*sqrt(4*x3^2*(c^2-a^2)+4*x3*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
y4=(1/(2*a))*sqrt(4*x4^2*(c^2-a^2)+4*x4*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
y5=-(1/(2*a))*sqrt(4*x1^2*(c^2-a^2)+4*x1*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2);
y6=-(1/(2*a))*sqrt(4*x2^2*(c^2-a^2)+4*x2*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
y7=-(1/(2*a))*sqrt(4*x3^2*(c^2-a^2)+4*x3*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
y8=-(1/(2*a))*sqrt(4*x4^2*(c^2-a^2)+4*x4*(a^2*c-c^3)+a^4+c^4-2*a^2*c^2); 
x5=x1;
x6=x2;
x7=x3;
x8=x4; 
x=[x1 x2 x3 x4 x5 x6 x7 x8]; 
y=[y1 y2 y3 y4 y5 y6 y7 y8]; 
 m=zeros(1,8);
 n=zeros(1,8); 
 for i=1:8
m(i)=-sqrt(x(i)^2+y(i)^2) + sqrt((x(i)-c)^2+y(i)^2);
 m(i)=real(m(i));    
 n(i)=-sqrt(x(i)^2+y(i)^2)+sqrt((x(i)-d)^2+(y(i)-e)^2);
  n(i)=real(n(i));
 end
at=(a-m(1))^2;
bt=(b-n(1))^2; 
dt = at + bt;   
xt = x(1); 
yt = y(1); 
for j=2:8
    ap = (a-m(j))^2;
    bp = (b-n(j))^2;
    de = ap + bp;    
    if(de<dt)
        xt=real(x(j));       
        yt=real(y(j));        
        dt = ap + bp;
    end
end 
res(1) = xt; %here is stored x*
res(2) = yt; %here is stored y*
end
