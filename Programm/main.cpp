#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

#define g 9.81 // N/kg

class Doppelpendel {
  public:
    double l1;
    double l2;
    double m1;
    double m2;

    double t;

    double phi1;
    double dphi1;
    double d2phi1;
    double phi2;
    double dphi2;
    double d2phi2;

    double x1;
    double y1;
    double x2;
    double y2;

    //Konstruktor
    Doppelpendel();

    //init
    void initalize(double, double, double, double, double, double);

    // Berechne die Trajektorie der beiden Körper
    void berechne_trajektorie();

    // Runge-Kutta 2.Ordnung
    void RK2_Integrator(double);

    // Runge-Kutta 4.Ordnung
    void RK4_Integrator(double);

    double d2phi1dt2 (double, double, double, double);
    double d2phi2dt2 (double, double, double, double);
};

Doppelpendel::Doppelpendel(){
    l1 = 1.0; //m
    m1 = 1.0; //kg
}

void Doppelpendel::initalize(double _l2, double _m2, double _phi1, double _dphi1dt, double _phi2, double _dphi2dt){
    t=0;

    l2 = _l2;
    m2 = _m2;
    phi1 = _phi1;
    phi2 = _phi2;
    dphi1 = _dphi1dt;
    dphi2 = _dphi2dt;
}


void Doppelpendel::berechne_trajektorie(){
    ofstream out;
    out.open("test.txt");
    double dt = 0.1;
    for(t=0; t<10000; t+=dt) {
        x1 = l1*sin(phi1);
        y1 = -l1*cos(phi1);
        x2 = x1+l2*sin(phi2);
        y2 = y1-l2*cos(phi2);
        out << t << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << endl;

        RK2_Integrator(dt);
    }

}

void Doppelpendel::RK2_Integrator(double dt) {
    double k1_1 = dt*dphi1;
    double k1_2 = dt*dphi2;
    double n1_1 = dt*d2phi1dt2(phi1,dphi1,phi2,dphi2);
    double n1_2 = dt*d2phi2dt2(phi1,dphi1,phi2,dphi2);

    double k2_1 = dt*(dphi1+n1_1/2.);
    double k2_2 = dt*(dphi1+n1_1/2.);
    double n2_1 = dt*d2phi1dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);
    double n2_2 = dt*d2phi2dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);

    phi1 += k2_1;
    phi2 += k2_2;
    dphi1 += n2_1;
    dphi2 += n2_2;
}

void Doppelpendel::RK4_Integrator(double dt) {

}

double Doppelpendel::d2phi1dt2(double _phi1, double _dphi1, double _phi2, double _dphi2){
    double c = cos(_phi1-_phi2);
    double s = sin(_phi1-_phi2);
    return -(m2*l1*_dphi1*_dphi1*s*c-m2*g*sin(_phi2)*c+m2*l2*_dphi2*_dphi2*s+g*(m1+m2)*sin(_phi1))/(l1*((m1+m2)-m2*c*c));
}

double Doppelpendel::d2phi2dt2(double _phi1, double _dphi1, double _phi2, double _dphi2){
    double c = cos(_phi1-_phi2);
    double s = sin(_phi1-_phi2);
    return ((m1+m2)*l1*_dphi1*_dphi1*s-(m1+m2)*g*sin(_phi2)+m2*l2*_dphi2*_dphi2*s*c+g*(m1+m2)*sin(_phi1))/(l2*((m1+m2)-m2*c*c));
}

int main()
{
    Doppelpendel pendel;
    double l2=10;
    double m2=0.1;
    pendel.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);
    pendel.berechne_trajektorie();
    return 0;
}

