%%PR Newton Raphson
tc=input('Critical temp in kelvin');
pc=input('Critical pressure in bar');
omega=input('Input omega');
t=input('Temp in kelvin');
p=input('Pressure in bar');
tr=t/tc;
a=0.45724*((0.08314*tc)^2)/pc;
b=0.07780*0.08314*tc/pc;
if omega<=0.49
    k=0.37464+1.54226*omega-0.26992*(omega^2);
else
    k=0.379642+1.48503*omega-0.164423*(omega^2)+0.016666*(omega^3);
end
alpha=1+k*(1-(tr^0.5)^2);
A=alpha*a*p/((0.08314*t)^2);
B=b*p/(0.08314*t);
%newton raphson
x0=1;
count=1;
array_x1= [];
array_x0=[];
m=0.001;
for i=1:100
    x1=((x0^3)-((1-B)*(x0^2))+(A-2*B-3*(B^2))*x0-((A*B)-(B^2)-(B^3)))/((3*x0^2)-((1-B)*(2*x0))+(A-2*B-3*(B^2)));
    array_x1=[array_x1,x1]; 

    x0=x0-x1;
    array_x0=[array_x0,x0];

    if abs(x1)<m
        fprintf('\nZ equals %f\n',x0);
        V=x0*0.08314*t/p/1000;
        fprintf('Vm = %f m3/gmol\n',V);
        if count==3
            break;
        else
            count=count+1;
            x0=x0/2;
        end
    else
        i=i+1;
    end
%disp (array_x1)
disp (array_x0)
end