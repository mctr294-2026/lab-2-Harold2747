#include "roots.hpp"

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    
    if ( ( (f(a) > 0) && (f(b) > 0) )  || ( (f(a) < 0) && (f(b) < 0)))
        {
            return false ;
        }
    
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

        midpoint = (a+b/2);
        fc = f(midpoint);
    }

    *root = midpoint;
    
    return true;
    // Ask if this method can ever fail DETERMINE LOCAL MAX
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    
    if ( ( (f(a)>0) && (f(b) > 0) )  || ( (f(a) < 0) && (f(b) < 0)))
        {
            return false ;
        }
    
    double midpoint = (f(a)*(b-a))/(f(b)-f(a));
    double fc = f(midpoint);
    
    while (fabs(0 - fc) > 1e-6)
    {
        if (f(a) * fc < 0)
        {
            b = midpoint;
        }

        if (f(b) * fc < 0)
        {
            a = midpoint;
        }

        midpoint = (f(a)*(b-a))/(f(b)-f(a));
        fc = f(midpoint);
    }
    
    *root = midpoint;

    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    
    // Is it okay to use multiple return values
    
    if (g(c) == 0)
    {
        return false;
    }

    double x1 = c;
    double x2 = c - f(c)/g(c);
    /*int iterations = 0;    Should I test a number of iterations or is just checking
    if it is within the range sufficient*/

    while ((fabs(x2-x1) > 1e-6))
    {
        if (g(x2) == 0 || (x2 < a) || (x2 > b))
        {
            return false;
        }

        double x3 = x2 - f(x2)/g(x2);
        x1 = x2;
        x2 = x3;
        // iterations++;
    }
    
   /*if (iterations == 100)
    {
        return false;
    } */
    
    *root = x2;

    return true;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    
    double SecondGuess = 2;  //WHAT IS SECOND GUESS VALUE
    double xn = SecondGuess;
    double x1 = c;
    
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

    *root = x2;
    
    return true;
}

