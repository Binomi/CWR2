#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>

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

    double vx1;
    double vy1;
    double vx2;
    double vy2;

    double E0;
    double E;

    double salto[2];

    //Konstruktor
    Doppelpendel();

    //init
    void initalize(double _l2, double _m2, double _phi1, double _dphi1dt, double _phi2, double _dphi2dt);

    void kartesisch();

    double energie();

    // Berechne die Trajektorie der beiden Körper
    void berechne_trajektorie(double t_max, double dt, bool useRK4);

    // Runge-Kutta 2.Ordnung
    void RK2_Integrator(double dt);

    // Runge-Kutta 4.Ordnung
    void RK4_Integrator(double dt);

    double d2phi1dt2(double _phi1, double _dphi1, double _phi2, double _dphi2);
    double d2phi2dt2(double _phi1, double _dphi1, double _phi2, double _dphi2);
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
    kartesisch();
    E0 = energie();
    salto[0] = 0;
    salto[1] = 0;
}

void Doppelpendel::kartesisch(){
    x1 = l1*sin(phi1);
    y1 = -l1*cos(phi1);
    x2 = x1+l2*sin(phi2);
    y2 = y1-l2*cos(phi2);

    vx1 = l1*cos(phi1)*dphi1;
    vy1 = l1*sin(phi1)*dphi1;
    vx2 = vx1+l2*cos(phi2)*dphi2;
    vy2 = vy1+l2*sin(phi2)*dphi2;
}

double Doppelpendel::energie(){
    return m1*g*y1+m2*g*y2+0.5*m1*(vx1*vx1+vy1*vy1)+0.5*m2*(vx2*vx2+vy2*vy2);
}


void Doppelpendel::berechne_trajektorie(double t_max, double dt, bool useRK4){
    /*
    ofstream out;
    stringstream l2str;
    l2str << l2;
    stringstream m2str;
    m2str << m2;
    stringstream t_max_str;
    t_max_str << t_max;
    stringstream dt_str;
    dt_str << dt;
    string datname = "sim_dat\\l2="+l2str.str()+"_m2="+m2str.str()+"_tmax="+t_max_str.str()+"_dt="+dt_str.str();
    if(useRK4)  datname+="_rk4.txt";
    else        datname+="_rk2.txt";
    out.open(datname.c_str());
    */
    t=0;
    while(t<t_max) {
        /*
        out << t << "  " << x1 << "  " << y1 << "  " << x2 << "  " << y2 << "  "
            << vx1 << "  " << vy1 << "  " << vx2 << "  " << vy2 << "  " << sqrt(vx2*vx2+vy2*vy2)  << "  " << E << "  "
            << phi1 << "  " << phi2 << endl;
        */
        double phi_alt[] = {phi1, phi2};
        if(useRK4) {
            RK4_Integrator(dt);
        } else {
            RK2_Integrator(dt);
        }
        kartesisch();
        E = energie();
        double phi_neu [] = {phi1, phi2};
        for(int i=0; i<2; i++) {
            if(floor(phi_neu[i]/(2.*M_PI)+0.5) != floor(phi_alt[i]/(2.*M_PI)+0.5))
                salto[i]++;
        }
    }
    //out.close();
    //CopyFileA(datname.c_str(), "test.txt",false);
}

void Doppelpendel::RK2_Integrator(double dt) {
    double k1_1 = dt*dphi1;
    double k1_2 = dt*dphi2;
    double n1_1 = dt*d2phi1dt2(phi1,dphi1,phi2,dphi2);
    double n1_2 = dt*d2phi2dt2(phi1,dphi1,phi2,dphi2);

    double k2_1 = dt*(dphi1+n1_1/2.);
    double k2_2 = dt*(dphi2+n1_2/2.);
    double n2_1 = dt*d2phi1dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);
    double n2_2 = dt*d2phi2dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);

    phi1 += k2_1;
    phi2 += k2_2;
    dphi1 += n2_1;
    dphi2 += n2_2;

    t+=dt;
}

void Doppelpendel::RK4_Integrator(double dt) {
    double k1_1 = dt*dphi1;
    double k1_2 = dt*dphi2;
    double n1_1 = dt*d2phi1dt2(phi1,dphi1,phi2,dphi2);
    double n1_2 = dt*d2phi2dt2(phi1,dphi1,phi2,dphi2);

    double k2_1 = dt*(dphi1+n1_1/2.);
    double k2_2 = dt*(dphi2+n1_2/2.);
    double n2_1 = dt*d2phi1dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);
    double n2_2 = dt*d2phi2dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);

    double k3_1 = dt*(dphi1+n2_1/2.);
    double k3_2 = dt*(dphi2+n2_2/2.);
    double n3_1 = dt*d2phi1dt2(phi1+k2_1/2.,dphi1+n2_1/2.,phi2+k2_2/2.,dphi2+n2_2/2.);
    double n3_2 = dt*d2phi2dt2(phi1+k2_1/2.,dphi1+n2_1/2.,phi2+k2_2/2.,dphi2+n2_2/2.);

    double k4_1 = dt*(dphi1+n3_1);
    double k4_2 = dt*(dphi2+n3_2);
    double n4_1 = dt*d2phi1dt2(phi1+k3_1,dphi1+n3_1,phi2+k3_2,dphi2+n3_2);
    double n4_2 = dt*d2phi2dt2(phi1+k3_1,dphi1+n3_1,phi2+k3_2,dphi2+n3_2);

    phi1 += (k1_1+2.*k2_1+2.*k3_1+k4_1)/6.;
    phi2 += (k1_2+2.*k2_2+2.*k3_2+k4_2)/6.;
    dphi1 +=(n1_1+2.*n2_1+2.*n3_1+n4_1)/6.;
    dphi2 +=(n1_2+2.*n2_2+2.*n3_2+n4_2)/6.;

    t+=dt;
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

//----------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------

double d1(Doppelpendel a, Doppelpendel b) {
    double dx1 = a.x1-b.x1;
    double dy1 = a.y1-b.y1;
    return sqrt(dx1*dx1+dy1*dy1);
}

double d2(Doppelpendel a, Doppelpendel b) {
    double dx2 = a.x2-b.x2;
    double dy2 = a.y2-b.y2;
    return sqrt(dx2*dx2+dy2*dy2);
}

int main()
{
    ofstream out("test24.txt");

    double l2=0.5;
    double m2=0.05;
    Doppelpendel pendel2, pendel4;
    pendel2.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);
    pendel4.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);

    double dt = 0.005;
    double d1sum = 0;
    double d2sum = 0;
    for(double t=0; t<20.; t+=dt) {
        pendel2.berechne_trajektorie(dt,dt,false);
        pendel4.berechne_trajektorie(dt,dt,true);
        d1sum += d1(pendel2, pendel4);
        d2sum += d2(pendel2, pendel4);
        out << t << "  " << d1sum << "  " << d2sum << "  "
                << pendel2.x1 << "  " << pendel2.y1 << "  " << pendel2.x2 << "  " << pendel2.y2 << "  "
                << pendel4.x1 << "  " << pendel4.y1 << "  " << pendel4.x2 << "  " << pendel4.y2 << "  "
                << pendel2.E << "  "  << pendel4.E  << endl;
    }
    out.close();
    cout << pendel2.salto[0] << "  " << pendel2.salto[1] << "\n" << pendel4.salto[0] << "  " << pendel4.salto[1] << endl;

    ofstream out2("testsplot.txt");
    for(int i=-10; i<10; i++) {
        for(int j=-10; j<10; j++) {
            out2 << i << "  " << j << "  " << i+j << endl;
            if(j%2!=0) out2 << endl;
        }
    }
    out2.close();
    return 0;
}
