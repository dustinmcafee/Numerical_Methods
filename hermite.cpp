#define _USE_MATH_DEFINES

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<stdio.h>
#include<iomanip>

using namespace std;

//double f(double x){return (100/pow(x, 2))*sin(10/x);}
//double df(double x){return (20 * pow(x, 3)) - (92.34 * pow(x, 2)) + (101.16 * x);}
double f(double x){return ((3*sin(x)) + (2*cos(-x)) + log(x) + x);}
double df(double x){return ((-2*sin(x)) + (3*cos(x)) + pow(x, -1) + 1);}

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

int main(int argc, char *argv[]){

	int n = stoi(argv[1]);

	vector<double> x (n);
	vector<double> y (n);
	vector<double> dy (n);
	vector<double> result;

	for(int i = 0; i < n; i++){
		x[i] = (double) (i+1);
		y[i] = f((double) (i+1));
		dy[i] = df((double) (i+1));
//	cout << x[i] << "  " << y[i] << "  " << dy[i] << endl;
	}

	result = hermite_interpolation(x, y, dy);

	cout << "The Coefficients for the Hermite Interpolation are: " << endl;
	for (int i = 0; i < result.size() - 2; i++){
		cout << i << " : " << result[i] << endl;
	}

	return 0;
}
