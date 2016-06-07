# Numerico-M-TODO-DE-BRENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                      %%%%%
%%%%%                       DEFINICIÓN DEL PROBLEMA:                       %%%%%
%%%%%    Mientras que el método de la secante y posición falsa convergen   %%%%%     
%%%%%      formalmente más rápido que el método de la bisección, uno       %%%%%
%%%%%     encuentra funciones patológicas de práctica para los cuales      %%%%%     
%%%%%   bisección converge más rápidamente. Estos pueden ser, funciones    %%%%%
%%%%%   discontinuas agitadas, o incluso funciones suaves si la segunda    %%%%%     
%%%%% derivada cambia abruptamente cerca de la raíz. Bisección siempre a   %%%%%
%%%%%  la mitad el intervalo, mientras que la secante y posición falsa a   %%%%%     
%%%%%    veces pueden pasar muchos ciclos acercándose lentamente a los     %%%%%
%%%%%    límites lejanos más cerca de la raíz. Por otro lado el Método     %%%%%     
%%%%%   Ridders hace un trabajo mucho mejor, pero también a veces puede    %%%%%
%%%%%       ser engañado. ¿Hay una manera de combinar la convergencia      %%%%%     
%%%%%            super-lineal con la seguridad de bisección?               %%%%%
%%%%%                                                                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                         VALORES DE ENTRADA:                     %%%%%
%%%%% func: Es una función a la cual se desea calcular una            %%%%%  
%%%%%      aproximación a su raíz.                                    %%%%%      
%%%%% x1,x2: Son los extremos del intervalo que contiene a la raíz    %%%%%
%%%%%      la función fun c                                           %%%%%   
%%%%% TOL: Es la tolerancia que deseamos tener en la aproximación     %%%%%   
%%%%%      de la raíz                                                 %%%%%   
%%%%%                                                                 %%%%%
%%%%%                         VALORES DE SALIDA:                      %%%%%
%%%%%             la función devuelve un vector [p itmax]             %%%%%     
%%%%% p: Es la aproximación a la raíz                                 %%%%%
%%%%% itmax: Es la cantidad de veces que se evaluo la función         %%%%%   
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p itmax] = Brent(func,x1,x2,TOL)


ITMAX=10000;        % valor máximo de iteraciones.
EPS=3.0*(10^-8);    % presición para punto flotante
a=x1;
b=x2;
fa=feval(func,a);  %Evalua la funcion en el punto dado
fb=feval(func,b);
if((fa>0) & (fb>0)) || ((fa<0) & (fb<0)); %Verifica si la raíz se encuentra en el intervalo dado.
    error('Procedimiento terminado sin exito, SELECIONE OTRO INTERVALO');
else
    c=b;
    fc=fb;
    iter=0;   % Guardar el número de iteraciones
    while iter<ITMAX;
        if((fb>0) & (fc>0)) || ((fb<0) & (fc<0));
            c=a; 
            fc=fa;
            d=(b-a);
            e=d;
        end
        if abs(fc)<abs(fb);
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        end
        tol1=2*EPS*abs(b)+ 0.5*TOL; %Verificación de Convergencia.
        xm=0.5*(c-b);
        if(abs(xm)<=tol1) || (fb==0);
            p=b;
            itmax=iter;
            break
        end
        if(abs(e)>=tol1 && abs(fa)>abs(fb))
            s=fb/fa; %Intentar interpolación cuadrática inversa.
            if(a==c);
                p=2*xm*s;
                q=1-s;
            else
                q=fa/fc;
                r=fb/fc;
                p=s*(2*xm*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            end
            if(p>0); %Comprobar si en los límites.
                q=-q;
                p=abs(p);
            end
            if(2*p < min(3*xm*q-abs(tol1*q),abs(e*q)));
                e=d; %Aceptar la interpolación.
                d=p/q;
            else
                d=xm; %La interpolación ha fallado. Usar el método de bisección.
                e=d;
            end
        else %Límites disminuyendo muy lentamente. Usar el método de bisección.
            d=xm;
            e=d;
        end
        a=b; %Mueva la última mejor estimación de a.
        fa=fb;
        if(abs(d)>tol1); %Evaluar la nueva raíz
            b=b+d;
        else
            if sign(xm)==0;
                b=b+abs(tol1)*(-1);
            else
                b=b+abs(tol1);
            end
        end
        fb=feval(func,b);
        iter=iter+1;
    end
    if iter==ITMAX
    	error ('La funcion Brent excedio el máximo de iteraciones');
    end
    p=b;
    itmax=iter;
end
