#include <iostream>
#include <math.h>
#include <iostream>

using namespace std;

#define G 9.81 // N/kg

class Doppelpendel {
  public:
    double l1;
    double l2;
    double m1;
    double m2;

    double phi1;
    double dphi1dt;
    double d2phi1dt2;
    double phi2;
    double dphi2dt;
    double d2phi2dt2;

    double x1;
    double y1;
    double x2;
    double y2;

    //Konstruktor
    Doppelpendel(double, double, double, double, double, double);

    // Berechne die Trajektorie der beiden Körper
    void berechne_trajektorie();

    // Runge-Kutta 2.Ordnung
    void RK2_Integrator(double);

    // Runge-Kutta 4.Ordnung
    void RK4_Integrator(double);
};

Doppelpendel::Doppelpendel(double _l2, double _m2, double _phi1, double _dphi1dt, double _phi2, double _dphi2dt){
    l1 = 1.0; //m
    m1 = 1.0; //kg

    l2 = _l2;
    m2 = _m2;
    phi1 = _phi1;
    phi2 = _phi2;
    dphi1dt = _dphi1dt;
    dphi2dt = _dphi2dt;
}

void Doppelpendel::berechne_trajektorie(){
    x1 = l1*sin(phi1);
    y1 = -l1*cos(phi1);
    x2 = x1+l2*sin(phi2);
    y2 = y1-l2*cos(phi2);
}


int main()
{
    cout << "Hello World!" << endl;
    return 0;
}

