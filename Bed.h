#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>

using std::cout;
using std::endl;

namespace pl = std::placeholders;

class Bed{
private:
    const double eps = 0.4;
    const double u = 0.01;
    const double a = 0.05;
    const double b = 0.3;
    const double qm = 0.0375;
    const double co = 1;
    const double rho_p = 1000;

    const double k_ads = 0.05;
    const double rho_const =  rho_p * (1 - eps) / eps;

    const int number_of_nodes = 101;
    const double length = 0.5;
    const double dz = length / (number_of_nodes - 1);
    const double Ci_feed = 1;


public:

    std::vector <double> z; // Axial direction co-cordinate
    std::vector <double> S; // Solution vector with the form: [C0, ..., Cn, q0, ..., qn]
    
    Bed(){
        z.resize(number_of_nodes);
        S.resize(2*number_of_nodes);

        for (size_t i = 0; i < z.size(); i++)
        {
            z[i] = i * dz;
            cout << "z = " << z[i] << "\n";
        }

        // Assign zeros for the initial condition
        for (size_t i = 0; i < S.size(); i++)
        {
            S[i] = 0.0;
        }
        S[0] = Ci_feed;

    }

    void adsorption(const std::vector<double> &s , std::vector<double> &dsdt , double t );
    //static void log_adsorption(const std::vector<double> &s , const double t );
    void run_simulation();

    static void log_adsorption(const std::vector<double> &s , const double t ){
        cout << "t: " << t << '\t' << "C: " << s[100] << '\t' << "q out: " << s[201] << endl;
    }

    /*
    void operator() (const std::vector<double> &s , std::vector<double> &dsdt , const double t )
    {
        // Inlet nodes
        dsdt[0 + number_of_nodes] = k_ads * ((qm * b * s[0])/(1 + b * s[0]) - s[0 + number_of_nodes]);
        dsdt[0] = Ci_feed; 

        // Middle and final nodes
        for (int i = 1; i < number_of_nodes; i++) 
        {
            dsdt[i + number_of_nodes] = k_ads * ((qm * b * s[i])/(1 + b * s[i]) - s[i + number_of_nodes]);
            dsdt[i] = -u / dz * (s[i] - s[i - 1]) - rho_const * dsdt[i + number_of_nodes];
        }
    }*/

    void operator() (const std::vector<double> &s , const double t ){
        cout << "t: " << t << '\t' << "C: " << s[number_of_nodes - 1] << '\t' << "q out: " << s[(2 * number_of_nodes) - 1] << endl;
    }

};


void Bed::adsorption(const std::vector<double> &s , std::vector<double> &dsdt , double t )
{
    // Inlet nodes
    dsdt[0 + number_of_nodes] = k_ads * ((qm * b * s[0])/(1 + b * s[0]) - s[0 + number_of_nodes]);
    dsdt[0] = 0; 

    // Middle and final nodes
    for (int i = 1; i < number_of_nodes; i++) 
    {
        dsdt[i + number_of_nodes] = k_ads * ((qm * b * s[i])/(1 + b * s[i]) - s[i + number_of_nodes]);
        dsdt[i] = -u / dz * (s[i] - s[i - 1]) - rho_const * dsdt[i + number_of_nodes];
    }
}

void Bed::run_simulation(){
    cout << "Simulation Starting..." << endl;
    //boost::numeric::odeint::integrate( std::bind(&Bed::adsorption, std::ref(*this), pl::_1, pl::_2, pl::_3), S , 0.0 , 200.0 , 0.01 , log_adsorption);
    boost::numeric::odeint::integrate( std::bind(&Bed::adsorption, std::ref(*this), pl::_1, pl::_2, pl::_3), S , 0.0 , 1000.0 , 0.001, log_adsorption);
    cout << "Simulation Finished!" << endl;
}
