#include "roots.hpp"

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    
    /* Saying that if both f(a) and f(b) are greater than or less than zero, then
    this method won't necessarily work as there may not be any point where the
    function crosses zer0 within the interval */
    if ( ( (f(a) > 0) && (f(b) > 0) )  || ( (f(a) < 0) && (f(b) < 0)))
        {
            return false ;
        }
    
    /* Performing the biscection method operation as many times as it takes until
    f(c) equals zero with an accuracy of 1e-6 */
    double midpoint = (a+b)/2 ;
    double fc = f(midpoint);

    while (fabs(0 - fc) > 1e-6)
    {
        if (fc < 0)
        {
            a = midpoint;
        }
        else 
        {
            b = midpoint;
        }

        midpoint = (a+b)/2;
        fc = f(midpoint);
    } 

    /* Defining the value at the address of root to be the calculated value of c */
    *root = midpoint;
    
    return true;
    
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    
    /* Saying that if both f(a) and f(b) are greater than or less than zero, then
    this method won't necessarily work as there may not be any point where the
    function crosses zer0 within the interval */
    if ( ( (f(a)>0) && (f(b) > 0) )  || ( (f(a) < 0) && (f(b) < 0)))
        {
            return false ;
        }
    
    /* Performing the regula falsi operation as many times as it takes until
    f(c) equals zero with an accuracy of 1e-6 */
    double midpoint = (a-(f(a)*(b-a))/(f(b)-f(a)));
    double fc = f(midpoint);
    
    int iterations = 0;
    while (fabs(0 - fc) > 1e-6 && (iterations<1e7))
    {
        if (f(a) * fc < 0)
        {
            b = midpoint;
        }

        if (f(b) * fc < 0)
        {
            a = midpoint;
        }

        midpoint = (a-(f(a)*(b-a))/(f(b)-f(a)));
        fc = f(midpoint);
        iterations++;
    }
    
    if (iterations == 1e7)
    {
        return false;
    }

    /* Defining the value at the address of root to be the calculated value of c */
    *root = midpoint;

    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {

    
    /* Saying that if the derivative of the function equals zero that this method will 
    not work */
    if (g(c) == 0)
    {
        return false;
    }

    /* Performing the Newton-Raphson's method as many times as it takes to have x2 - x1 equal
    0 with an accuracy of 1e-6. But if at any point the derivative of x2 equals zero or x2
    is less than a or greater than b, the operation will be stopped and a false value will be 
    returned */
    double x1 = c;
    double x2 = c - f(c)/g(c);
    int iterations = 0;

    while ((fabs(x2-x1) > 1e-6) && (iterations < 1e6))
    {
        if (g(x2) == 0 || (x2 < a) || (x2 > b))
        {
            return false;
        }

        double x3 = x2 - f(x2)/g(x2);
        x1 = x2;
        x2 = x3;
        iterations++;
    }
    
    if (iterations == 1e6)
    {
        return false;
    }
    
    /* Defining the value at the address of root to be the calculated value of x2 */
    *root = x2;

    return true;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    
    double SecondGuess = (a+b)/2;
    double xn = SecondGuess;
    double x1 = c;
    
    /* Performing the secant method operation as many times as it takes to get x2-x1 to
    equal zero with a certainty of 1e-6. If at any point x2 is greater than b or less than a
    the operation will be stopped and a false value will be returned*/
    double x2 = x1 - f(x1) * (x1 - xn)/(f(x1)-f(xn));

    while (fabs(x2-x1) > 0)
    {
        if (x2 > b || x2 <a)
        {
            return false;
        }

        double x3 = x2 - f(x2) * (x2 - x1)/(f(x2) - f(x1));

        x1 = x2;
        x2 = x3;
    }

    /* Defining the value at the address of root to be the calculated value of x2 */
    *root = x2;
    
    return true;
}

