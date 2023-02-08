#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// IT WOULD BE NICE TO SHOW RATE OF CONEVRGENCE - DO IT


double f(double x){
    return x*x*x*x;
}

double minus_f(double x){
    return -x*x*x*x;
}

double f_dash( double x){
    return 4*x*x*x;
}

double f_double_dash(double x){
    return 12*x*x;
}

double g(double x){
    return sqrt(1 + x*x);
}

double minus_g(double x){
    return -sqrt(1 + x*x);
}

double g_dash( double x){
    return x / sqrt(1 + x*x);
}

double g_double_dash(double x){
    return (sqrt(1 + x*x) - (x*x)/sqrt(1 + x*x)) / (1 + x*x);
}

//returns the x value for maxima.
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

//returns the maxima
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

//returns the minima
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
        e = (fabs(x - x_old) / (x+0.01) )*100;
        i++;
        printf("%d | %.15lf | %.15lf | %.15lf | %.15lf | %.15lf\n", i, x, f(x), f_dash(x), f_double_dash(x), e);
    }
    return x;
}

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
        else if (func(x1)<func(x2))
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

//gives maxima
double twobythree(double a, double b, double e, double (*func)(double x))
{
    double xl, xu,x1,x2,xopt;
    xl = a;
    xu = b;
    
    double R = 2.0/3.0;
    double d = R*(xu - xl);
    
    x1 = xl + d;
    
    x2 = xu - d;
   
    bool exite = false;
    
    double error;
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
    }
    return xopt; 
}

void main(){

    // for x^4
    
    printf("\n\nOptimisation of f(x) = -x^4 using golden search method:\n\n");
    double result0 = golden(-1, 2, 1, &minus_f);
    printf("x_max = %.15lf\ny_max = %.15lf\n", result0, f(result0));

    printf("\n\nOptimisation of f(x) = -x^4 using parabolic interpolation:\n\n");
    double result = parabolic_interpolation_sequential(&minus_f, -2,1,5,8);
    printf("x_max = %.15lf\ny_max = %.15lf\n", result, f(result));

    printf("\n\nOptimisation of f(x) = x^4 using newton's method [initial point = 0.9]:\n\n");
    double result2 = newton(&f, &f_dash, &f_double_dash, 0.9, 1);
    printf("x_min = %.15lf\ny_min = %.15lf\n", result2, f(result2));

    printf("\n\nOptimisation of f(x) = x^4 using newton's method [initial point = 1.1]:\n\n");
    double result3 = newton(&f, &f_dash, &f_double_dash, 1.1, 1);
    printf("x_min = %.15lf\ny_min = %.15lf\n", result3, f(result3));

    //for  sqrt(1 + x^2)

    // printf("\n\nOptimisation of f(x) = -sqrt(1+x^2) using golden search method:\n\n");
    // double result7 = golden(-1, 2, 1, &minus_g);
    // printf("x_max = %.15lf\ny_max = %.15lf\n", result7, f(result7));

    printf("\n\nOptimisation of g(x) = -sqrt(1+x^2) using parabolic interpolation:\n\n");
    double result4 = parabolic_interpolation_sequential(&minus_g, -1,1,4,8 );
    printf("x_max = %.15lf\ny_max = %.15lf\n", result4, g(result4));

    printf("\n\nOptimisation of g(x) = sqrt(1+x^2) using newton's method [initial point = 0.9]:\n\n");
    double result5 = newton(&g, &g_dash, &g_double_dash, 0.9, 1);
    printf("x_min = %.15lf\ny_min = %.15lf\n", result5, g(result5));

    printf("\n\nOptimisation of g(x) = sqrt(1+x^2) using newton's method [initial point = 1.1]:\n\n");
    double result6 = newton(&g, &g_dash, &g_double_dash, 1.1, 1);
    printf("x_min = %.15lf\ny_min = %.15lf\n", result6, g(result6));
}