#include <stdio.h>
#include <stdbool.h>
#include <math.h>

//gives maxima
double golden(double a, double b, double e, double (*func)(double x))
{
    double xl, xu,x1,x2,xopt;
    xl = a;
    xu = b;
    double R = (sqrt(5) - 1)/2;
    double d = R*(xu - xl);
    x1 = xl + d;
    x2 = xu - d;
    bool exite = false;
    int i = 0;
    double error;
    printf("i |       x_opt       |          e        \n"); 
    while(exite == false)
    {
        if (func(x1)>func(x2))
        {
            xopt = x1;
            error = (1-R)*fabs((xu-xl)/xopt);
            xl = x2;
            x2 = x1;
            d = R*(xu - xl);
            x1 = xl + d;
        }
        else
        {
            xopt = x2;
            error = (1-R)*fabs((xu-xl)/xopt);
            xu = x1;
            x1 = x2;
            d = R*(xu - xl);
            x2 = xu - d;
        }
           if(error<e)
           {
             exite = true;
           }
        i++;
        printf("%d | %.15lf | %.15lf \n", i, xopt, error);
    }
    return xopt;   
}

//return the maxima
double twobythree(double a, double b, double e, double (*func)(double x))
{
    double xl, xu,x1,x2,xopt;
    xl = a;
    xu = b;
    
    double R = 2.0/3.0;
    double d = R*(xu - xl);
    
    x1 = xl + d;
    
    x2 = xu - d;
    int i = 0;
    bool exite = false;
    
    double error;
    printf("i |       x_opt       |          e        \n"); 
    while(exite == false)
    {
        if (func(x1)>func(x2))
        {
            xopt = x1;
            error = (1-R)*fabs((xu-xl)/xopt);
            xl = x2;
           
            d = R*(xu - xl);
            x1 = xl + d;
            x2 = xu - d; 
        }
        else
        {
            xopt = x2;
            
            error = (1-R)*fabs((xu-xl)/xopt);
            xu = x1;
            
            d = R*(xu - xl);
            x1 = xl + d;
            x2 = xu - d;  
        }
           if(error<e)
           {
             exite = true;
           }
        i++;
        printf("%d | %.15lf | %.15lf \n", i, xopt, error);      
    }
    return xopt; 
}

//returns the maxima.
double parabolic_interpolation_bracketing_approach(double (*f)(double x), double x_0, double x_1, double x_2, int num_iter){
    double y_0, y_1, y_2, x_3, y_3;
    printf("i |          x_0      |       f(x_0)      |        x_1        |      f(x_1)       |        x_2        |      f(x_2)       |        x_3        |      f(x_3)       \n");
    for(int i = 0; i<num_iter; i++){
        y_0 = f(x_0);
        y_1 = f(x_1);
        y_2 = f(x_2);
        
        x_3 = ((y_0 * (x_1*x_1 - x_2*x_2)) + (y_1 * (x_2*x_2 - x_0*x_0)) + (y_2 * (x_0*x_0 - x_1*x_1))) / ((2*y_0*(x_1 - x_2)) + (2*y_1*(x_2 - x_0)) + (2*y_2*(x_0 - x_1)));
        y_3 = f(x_3);

        //for debugging
        printf("%d | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf\n", i+1, x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);

        //updating for next iter
        if(y_1 > y_3){
            if(x_1 < x_3){
                x_2 = x_3;
            } else {
                x_0 = x_3;
            }
        } else {
            if(x_1 < x_3){
                x_0 = x_1;
                x_1 = x_3;
            } else {
                x_2 = x_1;
                x_1 = x_3;
            }
        }
    }
    return x_3;
}

//return the maxima
double parabolic_interpolation_sequential(double (*f)(double x), double x_0, double x_1, double x_2, int num_iter){
    double y_0, y_1, y_2, x_3, y_3;
    printf("i |          x_0      |       f(x_0)      |        x_1        |      f(x_1)       |        x_2        |      f(x_2)       |        x_3        |      f(x_3)       \n");
    for(int i = 0; i<num_iter; i++){
        y_0 = f(x_0);
        y_1 = f(x_1);
        y_2 = f(x_2);
        
        x_3 = ((y_0 * (x_1*x_1 - x_2*x_2)) + (y_1 * (x_2*x_2 - x_0*x_0)) + (y_2 * (x_0*x_0 - x_1*x_1))) / ((2*y_0*(x_1 - x_2)) + (2*y_1*(x_2 - x_0)) + (2*y_2*(x_0 - x_1)));
        y_3 = f(x_3);

        //for debugging
        printf("%d | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf\n", i+1, x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);

        //updating for next iter
        x_0 = x_1;
        x_1 = x_2;
        x_2 = x_3;
    }
    return x_3;
}

//returns the extrema
double newton(double (*f)(double x), double (*f_dash)(double x), double (*f_double_dash)(double x), double x, double error_threshold){
    double x_old, e;
    int i = 1;
    x_old = x;
    e = 100;
    printf("i |          x        |        f(x)       |       f'(x)       |      f''(x)       |          e        \n");
    printf("%d | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf\n", i, x, f(x), f_dash(x), f_double_dash(x), e);
    while(fabs(e) > error_threshold){
        x_old = x;
        x -= f_dash(x)/f_double_dash(x);
        e = (fabs(x - x_old) / (x) )*100;
        i++;
        printf("%d | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf\n", i, x, f(x), f_dash(x), f_double_dash(x), e);
    }
    return x;
}
