#include <iostream>
#include <cmath>

using namespace std;

/// Class declaratoin
/// .H files (Header files)
class circle
{
    double d_;

    public:

        /// Constructor

        /// Declaration: this function exists.
        circle(double x);

        /// Member functions
        void area();
};

/// Definition: what this function does.
/// .C files (Source)
circle::circle(double x): d_(x)
{}

void circle::area()
{
    cout<< 3.14 * pow(d_, 2) / 4;
}

/// Application
/// Solver or utility
int main()
{
    cout<< "Hello!\n";

    circle c1(2);
    c1.area();
}
