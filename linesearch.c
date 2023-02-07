#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// IT WOULD BE NICE TO SHOW RATE OF CONEVRGENCE - DO IT


double f(double x){
    return x*x*x*x;
}

double f_dash( double x){
    return 4*x*x*x;
}

double f_double_dash(double x){
    return 12*x*x;
}

//returns the x value for maxima.
double parabolic_interpolation(double (*f)(double x), double x_0, double x_1, double x_2, int num_iter){
    double y_0, y_1, y_2, x_3, y_3;
    printf("Header - fill it\n");
    for(int i = 0; i<num_iter; i++){
        y_0 = f(x_0);
        y_1 = f(x_1);
        y_2 = f(x_2);
        
        x_3 = ((y_0 * (x_1*x_1 - x_2*x_2)) + (y_1 * (x_2*x_2 - x_0*x_0)) + (y_2 * (x_0*x_0 - x_1*x_1))) / ((2*y_0*(x_1 - x_2)) + (2*y_1*(x_2 - x_0)) + (2*y_2*(x_0 - x_1)));
        y_3 = f(x_3);

        //for debugging
        printf("%d | %lf | %lf | %lf | %lf | %lf | %lf | %lf | %lf\n", i+1, x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);

        //updating for next iter - check this
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

double newton(double (*f)(double x), double (*f_dash)(double x), double (*f_double_dash)(double x), double x, double error_threshold){
    double x_old, e;
    int i = 1;
    x_old = x;
    e = 100;
    printf("Header - fill it\n");
    printf("%d | %lf | %lf | %lf | %lf | %lf\n", i, x, f(x), f_dash(x), f_double_dash(x), e);
    while(e > error_threshold){
        x_old = x;
        x -= f_dash(x)/f_double_dash(x);
        e = (fabs(x - x_old) / x )*100;
        i++;
        printf("%d | %lf | %lf | %lf | %lf | %lf\n", i, x, f(x), f_dash(x), f_double_dash(x), e);
    }
    return x;
}

void main(){
    //double result = parabolic_interpolation(&f, 0,1,4,100);
    //printf("x_max = %lf\ny_max = %lf\n", result, f(result));

    double result2 = newton(&f, &f_dash, &f_double_dash, 0.9, 1);
    printf("x_max = %lf\ny_max = %lf\n", result2, f(result2));
}