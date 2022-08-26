#ifndef RIMODE_H_INCLUDED
#define RIMODE_H_INCLUDED

class Rimode
{
    private:
        double _pr;
    public:
        Rimode(double pr) : _pr(pr) {}

        void operator() (const state_type &x,state_type &dxdt,const double /* t */)
        {
            dxdt[0]=x[1];
            dxdt[1]=(gravity*sin(x[0])-(_pr-ball_radius*sin(x[0]))*x[3]*x[3]*cos(x[0]))/(1.4*ball_radius);
            dxdt[2]=x[3];
            dxdt[3]=(12./7.)*(ball_radius*x[1]*x[3]*cos(x[0]))/(_pr-ball_radius*sin(x[0]));
        }
};

Rimode rimodem(mpocket_r);
Rimode rimodec(cpocket_r);
boost::numeric::odeint::runge_kutta4<state_type> stepper;

#endif // RIMODE_H_INCLUDED
