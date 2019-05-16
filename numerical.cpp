#define _USE_MATH_DEFINES

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<stdio.h>
#include<iomanip>

using namespace std;

double f(double x){return (x - cos(x));}
double df(double x){return (1 + sin(x));}

double bisection_method(double p1, double p2, double tol, int& count, int& max_iteration)
{
   double x, tmp, a = p1, b = p2;
   x = 1.0;
   tmp = 0;	// To ensure the while statement starts.
                // prone to error: f(x) == 0.
   cout << endl;
   while (f(x) != 0.0 && (x - tmp) > tol && count <= max_iteration)
   {
      double tmpa = a;
      double tmpb = b;
      x = a + ((b - a)/2);
      tmp = a;
      if((f(x) * f(a)) > 0){a = x;}
      else{b = x;}
      count++;
      cout << "a = " << tmpa << " " << "b = " << tmpb << " ";
      cout << "p" << count << " = " << x << " ";
      cout << "error: " << x - tmp << endl;
   }
   return x;
}

double secant_method(double p1, double p2, double tol, int& count, int& max_iteration)
{
   double x, a = p1, b = p2;
   double fa, fb;

   cout << endl;
   while ((abs(b - a)) > tol && count <= max_iteration)
   {
      fa = f(a);
      fb = f(b);
      x = (b - ((fb * (b - a))/(fb - fa)));
      count++;

      cout << "a = " << a << " " << "b = " << b << " " << "p" << count << " = " << x;

      a = b;
      b = x;
      cout << " error: " << abs(b - a) << endl;
   }
   return x;
}

double mullers_method(double p0, double p1, double p2, double tol, int& count, int& max_iteration)
{
   double p, h1, h2, d1, d2, d, b, bigd, e1, e2, bige, h;
   h1 = p1 - p0;
   h2 = p2 - p1;
   d1 = (f(p1) - f(p0))/h1;
   d2 = (f(p2) - f(p1))/h2;
   d = ((d2 - d1)/(h2 + h1));

   cout << endl;
   while (count <= max_iteration)
   {
      b = d2 + (h2 * d);

      if (pow(b, 2.0) - (4 * f(p2) * d) < 0) { cerr << "bigd imaginary at p" << count << endl; return -1.0;}

      bigd = sqrt(pow(b, 2.0) - (4 * f(p2) * d));
      if (abs(b - bigd) < abs(b + bigd)) {bige = b + bigd;}
      else {bige = b - bigd;}
      h = (-2 * f(p2))/bige;
      p = p2 + h;

      cout << "p0: " << p0 << " p1: " << p1 << " p2: " << p2 << " p" << count << " = " << p << " error: " << abs(h) << endl;

      if (abs(h) < tol) {return p;}
      p0 = p1;
      p1 = p2;
      p2 = p;
      h1 = p1 - p0;
      h2 = p2 - p1;
      d1 = (f(p1) - f(p0))/h1;
      d2 = (f(p2) - f(p1))/h2;
      d = ((d2 - d1)/(h2 + h1));
      count++;
   }
   return -1.0;
}

double false_position(double p1, double p2, double tol, int& count, int& max_iteration)
{
   double x, tmp, a = p1, b = p2;
   count = 0;
   x = 1.0;
   tmp = 100.0;            // to ensure while loop starts
   double fa, fb, fx;

   cout << endl;
   while (abs(x - tmp) > tol && count <= max_iteration)
   {
      fa = f(a);
      fb = f(b);
      x = (b - (fb * (b - a)/(fb - fa)));
      fx = f(x);

      cout << "a = " <<	a << " b = " <<	b << " p" << count << "	= " << x;

      count++;
      if((fx * fb) < 0){a = b; fa = fb;}
      tmp = b;
      b = x;
      fb = fx;

     cout << " error: " << abs(x - tmp) << endl;
   }
   return x;
}

double newtons_method(double p1, double tol, int& count, int& max_iteration, bool iterations)
{
   double x = 100.0, fx, dfx, tmp = 0.0, a = p1;

   while(abs(x - tmp) > tol && count <= max_iteration)
   {
       cout << "p" << count << ": " << a;
       fx = f(a);
       dfx = df(a);
       x = (a - (fx/dfx));
       count++;
       tmp = a;
       a = x;
       cout << " error: " << abs(x - tmp) << endl;
   }
   return x;
}

vector<vector<double> > neville(vector<double> x_values, vector<double> y_values, double x)
{
   if (x_values.size() != y_values.size()) {cerr << " x.size !=	y.size ";}
   vector<vector<double> > nev (x_values.size() + 1, vector<double> (x_values.size() + 1));

   for (int i = 0; i < nev.size(); i++)
   {
      for (int j = 0; j < nev[i].size(); j++)
      {
         nev[i][j] = 0.0;
      }
   }

   for (int i = 0; i < x_values.size(); i++)
   {
      nev[i][0] = y_values[i];
   }

   for (int i = 1; i <= x_values.size(); i++)
   {
      for (int j = 1; j <= i; j++)
      {
         nev[i][j] = (((x - x_values[i-j]) * nev[i][j-1]) - ((x - x_values[i]) * nev[i-1][j-1]))/(x_values[i] - x_values[i-j]);
      }
   }
   return nev;
}

double lagrange(vector<double> x_values, vector<double> y_values, double x)
{
   double first, second, third, fourth, fifth, sixth;
   int order = x_values.size() - 1;
   if (order == 0) {return y_values[0];}
   else if (order == 1)
   {
      first = (y_values[0] * (x - x_values[1]))/(x_values[0] - x_values[1]);
      second = (y_values[1] * (x - x_values[0]))/(x_values[1] - x_values[0]);
      return first + second;
   }
   else if (order == 2)
   {
      first = (y_values[0] * (x	- x_values[1]) * (x - x_values[2])) /
         ((x_values[0] - x_values[1]) * (x_values[0] - x_values[2]));
      second = (y_values[1] * (x - x_values[0]) * (x - x_values[2])) /
         ((x_values[1] - x_values[0]) * (x_values[1] - x_values[2]));
      third = (y_values[2] * (x - x_values[0]) * (x - x_values[1])) /
         ((x_values[2] - x_values[0]) * (x_values[2] - x_values[1]));
      return first + second + third;
   }
   else if (order == 3)
   {
      first = (y_values[0] * (x - x_values[1]) * (x - x_values[2]) * (x - x_values[3])) /
         ((x_values[0] - x_values[1]) * (x_values[0] - x_values[2]) * (x_values[0] - x_values[3]));
      second = (y_values[1] * (x - x_values[0]) * (x - x_values[2]) * (x - x_values[3])) /
         ((x_values[1] - x_values[0]) * (x_values[1] - x_values[2]) * (x_values[1] - x_values[3]));
      third = (y_values[2] * (x	- x_values[0]) * (x - x_values[1]) * (x - x_values[3])) /
         ((x_values[2] - x_values[0]) * (x_values[2] - x_values[1]) * (x_values[2] - x_values[3]));
      fourth = (y_values[3] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2])) /
         ((x_values[3] - x_values[0]) * (x_values[3] - x_values[1]) * (x_values[3] - x_values[2]));
      return first + second + third + fourth;
   }
   else if(order == 4)
   {
      first = (y_values[0] * (x - x_values[1]) * (x - x_values[2]) * (x - x_values[3]) * (x - x_values[4])) /
         ((x_values[0] - x_values[1]) * (x_values[0] - x_values[2]) * (x_values[0] - x_values[3]) * (x_values[0] - x_values[4]));
      second = (y_values[1] * (x - x_values[0]) * (x - x_values[2]) * (x - x_values[3]) * (x - x_values[4])) /
         ((x_values[1] - x_values[0]) * (x_values[1] - x_values[2]) * (x_values[1] - x_values[3]) * (x_values[1] - x_values[4]));
      third = (y_values[2] * (x - x_values[0]) * (x - x_values[1]) * (x - x_values[3]) * (x - x_values[4])) /
         ((x_values[2] - x_values[0]) * (x_values[2] - x_values[1]) * (x_values[2] - x_values[3]) * (x_values[2] - x_values[4]));
      fourth = (y_values[3] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2]) * (x - x_values[4])) /
         ((x_values[3] - x_values[0]) * (x_values[3] - x_values[1]) * (x_values[3] - x_values[2]) * (x_values[3] - x_values[4]));
      fifth = (y_values[4] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2]) * (x - x_values[3])) /
         ((x_values[4] - x_values[0]) * (x_values[4] - x_values[1]) * (x_values[4] - x_values[2]) * (x_values[4] - x_values[3]));
      return first + second + third + fourth + fifth;
   }
   else if(order == 5)
   {
      first = (y_values[0] * (x - x_values[1]) * (x - x_values[2]) * (x - x_values[3]) * (x - x_values[4]) * (x - x_values[5])) /
         ((x_values[0] - x_values[1]) * (x_values[0] - x_values[2]) * (x_values[0] - x_values[3]) * (x_values[0] - x_values[4]) * (x_values[0] - x_values[5]));
      second = (y_values[1] * (x - x_values[0]) * (x - x_values[2]) * (x - x_values[3]) * (x - x_values[4]) * (x - x_values[5])) /
         ((x_values[1] - x_values[0]) * (x_values[1] - x_values[2]) * (x_values[1] - x_values[3]) * (x_values[1] - x_values[4]) * (x_values[1] - x_values[5]));
      third = (y_values[2] * (x - x_values[0]) * (x - x_values[1]) * (x - x_values[3]) * (x - x_values[4]) * (x - x_values[5])) /
         ((x_values[2] - x_values[0]) * (x_values[2] - x_values[1]) * (x_values[2] - x_values[3]) * (x_values[2] - x_values[4]) * (x_values[2] - x_values[5]));
      fourth = (y_values[3] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2]) * (x - x_values[4]) * (x - x_values[5])) /
         ((x_values[3] - x_values[0]) * (x_values[3] - x_values[1]) * (x_values[3] - x_values[2]) * (x_values[3] - x_values[4]) * (x_values[3] - x_values[5]));
      fifth = (y_values[4] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2]) * (x - x_values[3]) * (x - x_values[5])) /
         ((x_values[4] - x_values[0]) * (x_values[4] - x_values[1]) * (x_values[4] - x_values[2]) * (x_values[4] - x_values[3]) * (x_values[4] - x_values[5]));
      sixth = (y_values[5] * (x - x_values[0]) * (x - x_values[1]) * (x  - x_values[2]) * (x - x_values[3]) * (x - x_values[4])) /
         ((x_values[5] - x_values[0]) * (x_values[5] - x_values[1]) * (x_values[5] - x_values[2]) * (x_values[5] - x_values[3]) * (x_values[5] - x_values[4])); 
      return first + second + third + fourth + fifth + sixth;
   }

   cerr << "order too high, not enough code generated to solve problem." << endl;
   return -1;
}

vector<double> forward_difference(vector<double> x, vector<double> y)
{
   vector<vector<double> > table (x.size(), vector<double> (x.size()));
   for (int i = 0; i < table.size(); i++)
   {
      table[i][0] = y[i];
   }

   for (int i = 1; i < x.size(); i++)
   {
      for (int j = 1; j <= i; j++)
      {
         table[i][j] = (table[i][j-1] - table[i-1][j-1]) / (x[i] - x[i-j]);
      }
   }

   vector<double> diff (table.size());
   for (int i = 0; i < table.size(); i++)
   {
      diff[i] = table[i][i];
   }
   return diff;
}

vector<double> hermite_interpolation (vector<double> x, vector<double> y, vector<double> dy)
{
   int n = x.size() + 1;
   vector<double> z ((2*n)+2);
   vector<vector<double> > q ((2*n)+2, vector<double> ((2*n)+2));

   for (int i = 0; i <= n; i++)
   {
      z[2*i] = x[i];
      z[(2*i)+1] = x[i];
      q[2*i][0] = y[i];
      q[(2*i)+1][0] = y[i];
      q[(2*i)+1][1] = dy[i];
      
      if (i != 0)
      {
         q[2*i][1] = (q[2*i][0] - q[(2*i)-1][0]) / (z[2*i] - z[(2*i)-1]);
      }
   }
   for (int i = 2; i <= (2*n)+1; i++)
   {
      for (int j = 2; j <= i; j++)
      {
         q[i][j] = (q[i][j-1] - q[i-1][j-1]) / (z[i] - z[i-j]);
      }
   }

   vector<double> result (2*n);
   for (int i = 0; i < result.size(); i++)
   {
      result[i] = q[i][i];
   }
   return result;
}

double simpsonscomp(double a, double b, int n)
{
   double x = 0.0;
   double h = (b - a) / n;
   vector<double> xi (3);

   xi[0] = f(a) + f(b);
   xi[1] = 0;	// Summation of f(x[2*i-1]).
   xi[2] = 0;	// Summation of f(x[2*i]).

   for(int i = 1; i <= n - 1; i++)
   {
      x = a + (i * h);
      
      if(pow(-1, i) > 0)	// If i is even.
      {
         xi[2] += f(x);
      }
      else
      {
         xi[1] += f(x);
      }
   }

   double result = h * (xi[0] + (2 * xi[2]) + (4 * xi[1])) / 3;
   return result;
}

double trapezoidcomp(double a, double b, int n)
{
   double result;
   double h = (b - a) / n;
   double sum = 0.0;
   vector<double> x (n + 1);

   for (int i = 0; i < n; i++)
   {
      x[i] = a + (i * h);
   }
   x[n-1] = b;

   for (int i = 1; i <= n - 1; i++)
   {
      sum += f(x[i]);
   }
   result = (h / 2) * (f(a) + (2 * sum) + f(b));
   return result;
}

void romberg(double a, double b, int n)
{
   vector<vector<double> > r (n+10, vector<double> (n+10));
   double h = b - a;

   r[1][1] = (h/2)*(f(a)+f(b));

   cout << endl << fixed << setprecision(16) <<  r[1][1] << endl;

   for(int i = 2; i <= n; i++)
   {
      double sum = 0.0;

      for(int k = 1; k <= (pow(2, i - 2)); k++)
      {
         sum += f(a + ((k - 0.5) * h));
      }

      r[2][1] = 0.5 * (r[1][1] + (h * sum));					//Approximation from Trapezoidal Method.
      for(int j = 2; j <= i; j++)
      {
         r[2][j] = r[2][j-1] + ((r[2][j-1] - r[1][j-1]) / (pow(4, j-1) - 1));	//Extrapolation
      }

      h = h / 2;
      for(int j = 1; j <= i; j++)
      {
         cout << fixed << setprecision(16) << r[2][j] << "     ";
      }
      cout << endl;

      for(int j = 1; j <= i; j++)
      {
        r[1][j] = r[2][j];
      }
   }
} 

double adaptive_quadrature(double start, double end, double tolerance, int n)
{
	double approx = 0.0;
	double fd, fe, s1, s2;
	vector<double> v(9);
	int i = 1;

	vector<double> tol(n+10), a(n+10), h(n+10), fa(n+10), fc(n+10), fb(n+10), s(n+10), l(n+10);

	tol[i] = 10 * tolerance;
	a[i] = start;
	h[i] = (end - start) / 2;
	fa[i] = f(start);
	fc[i] = f(start + h[i]);
	fb[i] = f(end);
	s[i] = h[i] * (fa[i] + (4 * fc[i]) + fb[i]) / 3;		//Approximation from Simpson's method for entire interval.
	l[i] = 1;

	while (i > 0)
	{
		fd = f(a[i] + h[i] / 2);
		fe = f(a[i] + 3 * h[i] / 2);
		s1 = h[i] * (fa[i] + 4 * fd + fc[i]) / 6;			//Approximations from Simpson's method for halves of subintervals.
		s2 = h[i] * (fc[i] + 4 * fe + fb[i]) / 6;
		v[1] = a[i];										//Save data at this level;
		v[2] = fa[i];
		v[3] = fc[i];
		v[4] = fb[i];
		v[5] = h[i];
		v[6] = tol[i];
		v[7] = s[i];
		v[8] = l[i];

		i--;

		if (abs(s1 + s2 - v[7]) < v[6])
		{
			approx += s1 + s2;
		}
		else
		{
			if (v[8] >= n)
			{
				cout << "LEVEL EXCEEDED" << endl;			//Procedure Fails.
				return approx;
			}
			else											//Add one level.
			{
				i++;										//Data for right half of subinterval.
				a[i] = v[1] + v[5];
				fa[i] = v[3];
				fc[i] = fe;
				fb[i] = v[4];
				h[i] = v[5] / 2;
				tol[i] = v[6] / 2;
				s[i] = s2;
				l[i] = v[8] + 1;

				i++;										//Data for left half of subinterval.

				a[i] = v[1];
				fa[i] = v[2];
				fc[i] = fd;
				fb[i] = v[3];
				h[i] = h[i - 1];
				tol[i] = tol[i - 1];
				s[i] = s1;
				l[i] = l[i - 1];
			}
		}
	}
	return approx;											//Approximates the integral to within tolerance.
}









   
int main ()
{
   int test, max_iteration, count = 1;
   double p0, p1, tol, result;
   cout << "What test would you like to use?"
        << endl << "1. Bisection" << endl
        << "2. Secant Method" << endl
        << "3. Method of False Position" << endl
        << "4. Muller's Method" << endl 
        << "5. Newton's Method" << endl
        << "6. Neville's Method" << endl
        << "7. Lagrange Interpolation" << endl 
        << "8. Newton's Divided Difference Tables" << endl 
        << "9. Hermite Interpolation" << endl 
        << "10. Composite Simpson's Rule" << endl 
        << "11. Composite Trapezoidal Rule" << endl 
        << "12. Romberg" << endl 
		<< "13. Adaptive Quadrature" << endl << "Test = ";
   cin >> test;

   if (test != 6 && test != 7 && test != 8 && test != 9 && test != 10 && test != 11 && test != 12 && test != 13)
   {
      cout << "p0 = ";
      cin >> p0;
      if(test != 5)
      {
         cout << "p1 = ";
         cin >> p1;
      }
      cout << "tolerance = ";
      cin >> tol;
      cout << "Max Iterations = ";
      cin >> max_iteration;
   }

   if(test == 1)
   {
      if ((f(p0) * f(p1)) > 0) { cout << "a*b must be < 0" << endl; cin >> p0 >> p1; return -1;}
      result = bisection_method(p0, p1, tol, count, max_iteration);
   }
   else if(test == 2){result = secant_method(p0, p1, tol, count, max_iteration);}
   else if(test == 3){result = false_position(p0, p1, tol, count, max_iteration);}

   else if(test == 4)
   {
   	double p2;
   	cout << "p2 = ";
	cin >> p2;
	count = 3;
	result = mullers_method(p0, p1, p2, tol, count, max_iteration);
   }

   else if(test == 5)
   {
      cout << "Iterate for values from [a, b]? Y/N" << endl;
      char n;
      cin >> n;
      bool iterations = true;

      if (n == 'y' || n == 'Y')
      {
         iterations = false;
         double end, h, start = p0;
         cout << "a = p0" << endl;
         cout << "b = ";
         cin >> end;
         cout << endl << "Starting at a, how much to increment (h)?" << endl;
         cin >> h;

         for (double i = start; i <= end; i += h)
         {
            count = 1;
            cout << endl << "initial guess (p0) = " << i << endl;
            result = newtons_method(i, tol, count, max_iteration, iterations);

            if (count >= max_iteration) {cout << "The result most likely diverged. After "
                                              << count << " iterations, the result is p" 
                                              << count << " = " << result << endl;}
            else {cout << "The result is p" << count << " = " << result << endl;}
         }
      }
      else
      {
      result = newtons_method(p0, tol, count, max_iteration, iterations);
      }
   }

   else if(test == 6)    // Print Table for Lagrange Interpolation 
   {
      vector<double> x_values, y_values;      
      double input;
      const int SENTINEL = -1;

      cout << "x values to interpolate:" << endl;
      cout << "END WITH -1" << endl;
      while (input != SENTINEL)
      {
         cin >> input;
         x_values.push_back(input);
      }

      input = 0;
      cout << "Associated f(x) values: " << endl;
      cout << "END WITH -1" << endl;
      while (input != SENTINEL)
      {
         cin >> input;
         y_values.push_back(input);
      }

      cout << "Value to Approximate by Interpolation? x = ";
      input = 0;
      cin >> input;

      vector<vector<double> > nev = neville(x_values, y_values, input);

      for (int i = 0; i < x_values.size() - 1; i++)
      {
         cout << x_values[i] << " ";
         for (int j = 0; j <= i; j++)
         {
            cout << "P = " << nev[i][j] << " ";
         }
         cout << endl;
      }
   return 0;
   } 

   else if(test == 7)
   {
      vector<double> x_values, y_values;
      int order;
      double input;
      cout << "Order of Polynomial: " << endl;
      cin >> order;
      for(int i = 0; i <= order; i++)
      {
         cout << "x" << i << " = " << endl;
         cin >> input;
         x_values.push_back(input);
         cout << "f(x" << i << ") = " << endl;
         cin >> input;
         y_values.push_back(input);
      }
      cout << "X Value To Approximate = " << endl;
      cin >> input;
      result = lagrange(x_values, y_values, input);
   }

   else if (test == 8)
   {
      cout << "number of x values: " << endl;
      int n; cin >> n;
      vector<double> x (n);
      cout << "x values: " << endl;
      for (int i = 0; i < n; i++)
      {
         cin >> x[i];
      }

      cout << "Corresponding y values: " << endl;
      vector<double> y (n);
      for (int i = 0; i < n; i++)
      {
         cin >> y[i];
      }

      vector<double> diff = forward_difference(x, y);
      for (int i = 0; i < diff.size(); i++)
      {
         cout << i << " divided difference: " << diff[i] << endl;
      }
      return 0;
   }

   else if (test == 9)
   {
      cout << "Number of X values? (n + 1 = ?): " << endl;
      int z;
      cin >> z;

      vector<double> x (z);
      vector<double> y (z);
      vector<double> dy (z);
      vector<double> result;

      cout << "X values: " << endl;

      for (int i = 0; i < z; i++)
      {
         cin >> x[i];
      }

      cout << "Corresponding Y values: " << endl;

      for (int i = 0; i < z; i++)
      {
         cin >> y[i];
      }
      

      cout << "Corresponding DY/DX values: " << endl;

      for (int i = 0; i < z; i++)
      {
         cin >> dy[i];
      }

      result = hermite_interpolation(x, y, dy);

      cout << "The Coefficients for the Hermite Interpolation are: " << endl;

      for (int i = 0; i < result.size(); i++)
      {
         cout << i << " : " << result[i] << endl;
      }
      return 0;
   }
   
   else if (test == 10)
   {
      double a, b, result;
      int n;

      cout << "a = ";
      cin >> a;
      cout << "b = ";
      cin >> b;
      cout << "n = ";
      cin >> n;
      cout << endl;
      
      result = simpsonscomp(a, b, n);

      cout << "Integral of f(x) = " << setprecision(10) << result << endl;
  
      return 0;
   }

   else if (test == 11)
   {
      double a, b, result;
      int n;

      cout << "a = ";
      cin >> a;
      cout << "b = ";
      cin >> b;
      cout << "n = ";
      cin >> n;
      cout << endl;

      result = trapezoidcomp(a, b, n);

      cout << "Integral of f(x) = " << setprecision(10) << result << endl;

      return 0;
   }
   else if(test == 12)
   {
      double a, b;
      int n;

      cout << "a = ";
      cin >> a;
      cout << "b = ";
      cin >> b;
      cout << "n = ";
      cin >> n;
      romberg(a, b, n);
      return 0;
   }
   else if (test == 13)
   {
	   double a, b, tol, approx;
	   int n;

	   cout << "a = ";
	   cin >> a;
	   cout << "b = ";
	   cin >> b;
	   cout << "tolerance = ";
	   cin >> tol;
	   cout << "max number of levels = ";
	   cin >> n;

	   approx = adaptive_quadrature(a, b, tol, n);

	   cout << "approximation = " << setprecision(10) << approx << endl;

	   return 0;
   }


   else {cerr << "fail" << endl; return -1;}
   if ( count >= max_iteration && test != 6 && test !=7)
   {
      cout << "function most likely diverged/tolerance set too high Max Iterations == Count" << endl;
      return -1;
   }
   cout << "The result is p" << count << " = " << result << endl;
   return 0;
}
