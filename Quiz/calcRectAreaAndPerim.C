#include <iostream>
#include <cmath>

using namespace std;

/// Class declaratoin
/// .H files (Header files)
class rect
{
    double l_;
    double w_;

    public:

        /// Constructor

        /// Declaration: this function exists.
        rect(double, double);

        /// Member functions
        void area();
        void perim();
};

/// Definition: what this function does.
/// .C files (Source)
rect::rect(double a, double b): l_(a), w_(b)
{}

void rect::area()
{
    cout<< l_ * w_
        << "\n";
}

void rect::perim()
{
    cout<< 2 * (l_ + w_)
        << "\n";
}

/// Application
/// Solver or utility
int main()
{
    rect r(2,3); // Area = 6, P = 10

    cout<< "Area = ";
    r.area();

    cout<< "Perimeter = ";
    r.perim();
}
