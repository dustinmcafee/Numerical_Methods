#include<iostream>
#include<cmath>
#include<iomanip>

using namespace std;

double f(double x){return (x - cos(x));}
double df(double x){return (1 + sin(x));}

int main(int argc, char *argv[])
{

    int i=0,j=0,k;
    int n = stoi(argv[1]);
    double x[n], y[n][n];

    for (i=0;i<n;i++){            //input x values
        x[i] = (double) (i+1);
    }
    for (i=0;i<n;i++){
        y[i][0] = f((double) (i+1));
    }
    for (j=1;j<n;j++){            //loop to calculate the difference and store them in the matrix
        for (i=0;i<n-j;i++){
            y[i][j]=y[i+1][j-1]-y[i][j-1];
        }
    }

    cout.precision(2);        //set precision
    cout.setf(ios::fixed);
    cout<<"\n The Forward Difference Table is as follows: \n\n";
    cout<<"x"<<setw(10)<<"y"<<setw(10);    //formatting the output and creating table headings

    for (i=1;i<n;i++){
        cout<<"y"<<i<<setw(10);
    }

    cout<<"\n-----------------------------------------------------------------------\n";

    k=n;
    for (i=0;i<n;i++){            //loop for printing the diagonal matrix on the screen
        cout<<x[i]<<setw(10);
        for (j=0;j<k;j++){
            cout<<y[i][j];
            cout<<setw(10);
        }
        cout<<"\n";
        k--;
    }
return 0;
}
