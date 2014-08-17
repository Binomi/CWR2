#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;


// Erdbeschleunigung
#define g 9.81 // N/kg

class Doppelpendel {
public:
    //Längen und Massen
    double l1;
    double l2;
    double m1;
    double m2;

    //Zeit
    double t;

    //Winkel und Winkelgeschwindigkeiten
    double phi1;
    double dphi1;
    double phi2;
    double dphi2;

    //Kartesische Koordinaten
    double x1;
    double y1;
    double x2;
    double y2;
    //Geschwindigkeiten
    double vx1;
    double vy1;
    double vx2;
    double vy2;

    //Anfangsenergie und aktuelle Gesamtenergie
    double E0;
    double E;

    //Zähler der Überschläge
    double salto[2];

    //Konstruktor
    Doppelpendel();

    //Funktionen
    // init
    void initalize(double _l2, double _m2, double _phi1, double _dphi1dt, double _phi2, double _dphi2dt);

    // aus den Winkeln werden die kartesischen Koordinaten samt Ableitungen berechnet
    void kartesisch();

    // berechnet die aktuelle Gesamtenergie des Systems
    double energie();

    // Berechne die Trajektorie der beiden Körper
    void berechne_next_step(double dt, bool useRK4);

    // Runge-Kutta 2.Ordnung
    void RK2_Integrator(double dt);

    // Runge-Kutta 4.Ordnung
    void RK4_Integrator(double dt);

    // geben die zweiten Ableitungen zurück in Abh. von den Auslenkungen und Winkelgeschwindigkeiten
    double d2phi1dt2(double _phi1, double _dphi1, double _phi2, double _dphi2);
    double d2phi2dt2(double _phi1, double _dphi1, double _phi2, double _dphi2);
};

Doppelpendel::Doppelpendel(){
    l1 = 1.0; //m
    m1 = 1.0; //kg
}

void Doppelpendel::initalize(double _l2, double _m2, double _phi1, double _dphi1dt, double _phi2, double _dphi2dt){
    //reset
    t=0;
    salto[0] = 0;
    salto[1] = 0;

    //zweite Länge und Masse
    l2 = _l2;
    m2 = _m2;

    //Anfangsbedingungen
    phi1 = _phi1;
    phi2 = _phi2;
    dphi1 = _dphi1dt;
    dphi2 = _dphi2dt;

    kartesisch();

    //setzen der Anfangsenergie
    E0 = energie();
    E = E0;
}

void Doppelpendel::kartesisch(){
    //Koordinaten
    x1 = l1*sin(phi1);
    y1 = -l1*cos(phi1);
    x2 = x1+l2*sin(phi2);
    y2 = y1-l2*cos(phi2);

    //Geschwindigkeiten
    vx1 = l1*cos(phi1)*dphi1;
    vy1 = l1*sin(phi1)*dphi1;
    vx2 = vx1+l2*cos(phi2)*dphi2;
    vy2 = vy1+l2*sin(phi2)*dphi2;
}

double Doppelpendel::energie(){
    // Summe aller potentiellen und kinetischen Energien
    return m1*g*y1+m2*g*y2+0.5*m1*(vx1*vx1+vy1*vy1)+0.5*m2*(vx2*vx2+vy2*vy2);
}


void Doppelpendel::berechne_next_step(double dt, bool useRK4){

    //alte Winkelstellung
    double phi_alt[] = {phi1, phi2};

    //Integrationsschritt der Schrittweite dt
    if(useRK4) {
        RK4_Integrator(dt);
    } else {
        RK2_Integrator(dt);
    }

    // neue kartesische Koordinaten
    kartesisch();
    // neue Energie
    E = energie();

    //neue Winkelstellung
    double phi_neu [] = {phi1, phi2};
    // falls ...-pi,pi,3pi,... überschritten wird,
    // also der der jeweils höchste Punkt, wird der zugehörige Zähler erhöht
    for(int i=0; i<2; i++) {
        if(floor(phi_neu[i]/(2.*M_PI)+0.5) != floor(phi_alt[i]/(2.*M_PI)+0.5))
            salto[i]++;
    }
}

void Doppelpendel::RK2_Integrator(double dt) {
    //1.Zwischenschritt
    double k1_1 = dt*dphi1;
    double k1_2 = dt*dphi2;
    double n1_1 = dt*d2phi1dt2(phi1,dphi1,phi2,dphi2);
    double n1_2 = dt*d2phi2dt2(phi1,dphi1,phi2,dphi2);

    //2.Zwischenschritt
    double k2_1 = dt*(dphi1+n1_1/2.);
    double k2_2 = dt*(dphi2+n1_2/2.);
    double n2_1 = dt*d2phi1dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);
    double n2_2 = dt*d2phi2dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);

    //Werte aktualisieren
    t+=dt;
    phi1 += k2_1;
    phi2 += k2_2;
    dphi1 += n2_1;
    dphi2 += n2_2;
}

void Doppelpendel::RK4_Integrator(double dt) {
    //1.Zwischenschritt
    double k1_1 = dt*dphi1;
    double k1_2 = dt*dphi2;
    double n1_1 = dt*d2phi1dt2(phi1,dphi1,phi2,dphi2);
    double n1_2 = dt*d2phi2dt2(phi1,dphi1,phi2,dphi2);

    //2.Zwischenschritt
    double k2_1 = dt*(dphi1+n1_1/2.);
    double k2_2 = dt*(dphi2+n1_2/2.);
    double n2_1 = dt*d2phi1dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);
    double n2_2 = dt*d2phi2dt2(phi1+k1_1/2.,dphi1+n1_1/2.,phi2+k1_2/2.,dphi2+n1_2/2.);

    //3.Zwischenschritt
    double k3_1 = dt*(dphi1+n2_1/2.);
    double k3_2 = dt*(dphi2+n2_2/2.);
    double n3_1 = dt*d2phi1dt2(phi1+k2_1/2.,dphi1+n2_1/2.,phi2+k2_2/2.,dphi2+n2_2/2.);
    double n3_2 = dt*d2phi2dt2(phi1+k2_1/2.,dphi1+n2_1/2.,phi2+k2_2/2.,dphi2+n2_2/2.);

    //4.Zwischenschritt
    double k4_1 = dt*(dphi1+n3_1);
    double k4_2 = dt*(dphi2+n3_2);
    double n4_1 = dt*d2phi1dt2(phi1+k3_1,dphi1+n3_1,phi2+k3_2,dphi2+n3_2);
    double n4_2 = dt*d2phi2dt2(phi1+k3_1,dphi1+n3_1,phi2+k3_2,dphi2+n3_2);

    //Werte aktualisieren
    t+=dt;
    phi1 += (k1_1+2.*k2_1+2.*k3_1+k4_1)/6.;
    phi2 += (k1_2+2.*k2_2+2.*k3_2+k4_2)/6.;
    dphi1 +=(n1_1+2.*n2_1+2.*n3_1+n4_1)/6.;
    dphi2 +=(n1_2+2.*n2_2+2.*n3_2+n4_2)/6.;
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

// Abstand zwischen den ersten Massen zweier Doppelpendel
double d1(Doppelpendel a, Doppelpendel b) {
    double dx1 = a.x1-b.x1;
    double dy1 = a.y1-b.y1;
    return sqrt(dx1*dx1+dy1*dy1);
}

// Abstand zwischen den zweiten Massen zweier Doppelpendel
double d2(Doppelpendel a, Doppelpendel b) {
    double dx2 = a.x2-b.x2;
    double dy2 = a.y2-b.y2;
    return sqrt(dx2*dx2+dy2*dy2);
}

void variiereVerhaeltnisseQuantitativ(double l2_min, double l2_max, double l2_delta,
                                      double m2_min, double m2_max, double m2_delta,
                                      double T, double dt, string datname) {
    ofstream out(datname.c_str());
    Doppelpendel pendel2, pendel4;
    for(double l2=l2_min; l2<=l2_max; l2+=l2_delta){
        for(double m2=m2_min; m2<=m2_max; m2+=m2_delta) {
            pendel2.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);
            pendel4.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);
            double d1sum = 0.;
            double d2sum = 0.;
            for(double t=0; t<T; t+=dt) {
                pendel2.berechne_next_step(dt,false);
                pendel4.berechne_next_step(dt,true);
                d1sum += d1(pendel2, pendel4);
                d2sum += d2(pendel2, pendel4);
            }
            out << l2 << "  " << m2 << "  "
                << pendel2.salto[0] << "  " << pendel2.salto[1] << "  " << pendel4.salto[0] << "  " << pendel4.salto[1] << "  "
                << d1sum << "  " << d2sum << "  " << pendel2.E0 << "  " << pendel2.E << "  " << pendel4.E << endl;
        }
        out << endl;
    }
    out.close();
}

void rk2VSrk4(double l2, double m2, double T, double dt, string datname) {
    ofstream out(datname.c_str());
    out << "# " << l2 << "  " << m2 << endl;

    Doppelpendel pendel2, pendel4;
    pendel2.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);
    pendel4.initalize(l2,m2,M_PI/4.,1.0,M_PI/3.,1.5);

    double d1sum = 0.;
    double d2sum = 0.;

    for(double t=0; t<=T; t+=dt) {
        out << t << "  "
            << pendel2.x1 << "  " << pendel2.y1 << "  " << pendel2.x2 << "  " << pendel2.y2 << "  "
            << pendel2.vx1 << "  " << pendel2.vy1 << "  " << pendel2.vx2 << "  " << pendel2.vy2 << "  "
            << pendel4.x1 << "  " << pendel4.y1 << "  " << pendel4.x2 << "  " << pendel4.y2 << "  "
            << pendel4.vx1 << "  " << pendel4.vy1 << "  " << pendel4.vx2 << "  " << pendel4.vy2 << "  "
            << pendel2.E << "  "  << pendel4.E  << "  " << d1sum << "  " << d2sum << "  "
            << pendel2.salto[0] << "  " << pendel2.salto[1] << "  " << pendel4.salto[0] << "  " << pendel4.salto[1] << endl;

        pendel2.berechne_next_step(dt,false);
        pendel4.berechne_next_step(dt,true);
        d1sum += d1(pendel2, pendel4);
        d2sum += d2(pendel2, pendel4);
    }
    out.close();
}

int main()
{
    variiereVerhaeltnisseQuantitativ(0.01, 10., 0.01, 0.01, 10., 0.01, 10, 0.01, "l2m2_genauer1.txt");
    variiereVerhaeltnisseQuantitativ(0.01, 10., 0.01, 0.01, 10., 0.01, 20, 0.01, "l2m2_genauer2.txt");
    variiereVerhaeltnisseQuantitativ(0.01, 10., 0.01, 0.01, 10., 0.01, 30, 0.01, "l2m2_genauer3.txt");
    variiereVerhaeltnisseQuantitativ(0.01, 10., 0.01, 0.01, 10., 0.01, 40, 0.005, "l2m2_genauer4.txt");

    /*rk2VSrk4(5.75,7.,20.,0.005,"rk2rk4_1.txt");
    rk2VSrk4(0.5,0.05,20.,0.005,"rk2rk4_2.txt");
    rk2VSrk4(1.,1.,20.,0.005,"rk2rk4_3.txt");
    rk2VSrk4(3.,0.6,20.,0.005,"rk2rk4_4.txt");
    rk2VSrk4(3.,0.6,40.,0.005,"rk2rk4_5.txt");*/
    return 0;
}
