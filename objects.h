#ifndef OBJECTS_H_INCLUDED
#define OBJECTS_H_INCLUDED

#include <math.h>
#include <Eigen/Dense>
#include <SFML/Graphics.hpp>
#include <set>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include "polystuff.h"

using namespace boost::numeric::odeint;

//lengths are in inches.

const int framerate=100;
const double timestep=1./framerate;

const double epsilon=pow(10,-10);

const int updates=floor(1/(timestep*framerate));

const double in2m=1/39.37;
const double pi=3.1415926535897932384626433832795028841971;
const double gravity=9.81/in2m;

const double table_length=140.25;
const double table_width=70.125;
const double cush_thickness=1.875;
const double rail_thickness=2.445;
const double bounciness=0.7;
const double mus=0.0075; //0.005-0.015
const double muk=0.15; //0.15-0.4

const double dfactor=10;
const double window_width=(2*rail_thickness+2*cush_thickness+table_length)*dfactor;
const double window_height=(2*rail_thickness+2*cush_thickness+table_width)*dfactor;

const double ball_radius=0.02625/in2m;
const double ball_mass=0.141;
const double eball=1; //coefficient of restitution.

const double mpockets[2][2]={{rail_thickness+cush_thickness+table_width,rail_thickness+2*cush_thickness+table_width+0.156},
                             {rail_thickness+cush_thickness+table_width,rail_thickness-0.156}};
const double cpockets[4][2]={{rail_thickness,rail_thickness},
                             {rail_thickness,rail_thickness+2*cush_thickness+table_width},
                             {rail_thickness+2*cush_thickness+table_length,rail_thickness},
                             {rail_thickness+2*cush_thickness+table_length,rail_thickness+2*cush_thickness+table_width}};

const double spot_r=0.25;

const double mpocket_r=2.094;
const double cpocket_r=3.5;

const double pround=3.1;
const double k1=8.595-9.*sin(pi/4.);
const double k2=k1+4.5*sin(pi/4.)-2.445;

const double cueball_break_x=rail_thickness+cush_thickness+table_length-29;
const double cueball_break_y=rail_thickness+cush_thickness+table_width/2+6;

const double yellow_x=rail_thickness+cush_thickness+table_length-29;
const double yellow_y=rail_thickness+cush_thickness+table_width/2+11.687;
const double green_x=yellow_x;
const double green_y=rail_thickness+cush_thickness+table_width/2-11.687;
const double brown_x=yellow_x;
const double brown_y=rail_thickness+cush_thickness+table_width/2;
const double blue_x=rail_thickness+cush_thickness+table_width;
const double blue_y=brown_y;
const double pink_x=rail_thickness+cush_thickness+table_width/2;
const double pink_y=brown_y;
const double black_x=rail_thickness+cush_thickness+12.75;
const double black_y=brown_y;

//colours.
const int railcolour[3]={55,18,0};
const int baizecolour[3]={0,110,0};
const int cushioncolour[3]={0,80,0};
const int leathercolour[3]={255,229,153};

//grid squares.
std::array<std::array<std::array<int,7>,39>,73> grid={};
std::array<std::array<int,39>,73> grid_index={};

//cushion height-radius.
const double cushion_diff=0.03/in2m-ball_radius;
const double cushion_alpha=asin(cushion_diff/ball_radius);

//pocket things.
const double cpocket_angle=atan2(3.33402,1.06503);
const double mpocket_angle=atan2(1.93152,0.652752+0.156);

std::array<std::array<double,5>,3> cp={{{0.5*k1,0.5*k1,sqrt(2.)*0.5*k1-ball_radius,0.75*pi,-0.75*pi},
                                        {6.15,-2.445,4.5+ball_radius,-0.25*pi,-atan(6./11.)},
                                        {5.43,-1.125,3.+ball_radius,-atan(6./11.),0}}};
std::array<std::array<double,5>,3> mp={{{0,1.875,pround-ball_radius,pi-asin(1.625/pround),pi},
                                        {3.719,-0.438,2.094+ball_radius,-0.5*pi,-atan(0.625/0.812)},
                                        {4.344,-1.25,3.125+ball_radius,-atan(0.625/0.812),0}}};

std::array<double,2> xlim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_length-ball_radius};
std::array<double,2> ylim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_width-ball_radius};

std::array<double,2> mpline={blue_x-1.625,blue_x+1.625};
std::array<std::array<double,4>,8> cpline={{{1.,-k1,rail_thickness+k1,rail_thickness+k2},
                                            {1.,k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2*rail_thickness+2*cush_thickness+table_width-k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2*rail_thickness+2*cush_thickness+table_width+k1,rail_thickness+k1,rail_thickness+k2},
                                            {-1.,2*rail_thickness+2*cush_thickness+k1+table_length,rail_thickness+2*cush_thickness+table_length-k2+k1,rail_thickness+2*cush_thickness+table_length},
                                            {-1.,2*rail_thickness+2*cush_thickness-k1+table_length,rail_thickness+2*cush_thickness+table_length-k2,rail_thickness+2*cush_thickness+table_length-k1},
                                            {1.,-k1-table_width,rail_thickness+2*cush_thickness+table_length-k2+k1,rail_thickness+2*cush_thickness+table_length},
                                            {1.,k1-table_width,rail_thickness+2*cush_thickness+table_length-k2,rail_thickness+2*cush_thickness+table_length-k1}}};

//matrix.
Eigen::Matrix<double,46,46> M_;

std::tuple<double,double> add_vectors(double v1,double a1,double v2,double a2)
{
    double x=v1*sin(a1)+v2*sin(a2);
    double y=v1*cos(a1)+v2*cos(a2);

    return std::make_tuple(nroot(pow(x,2)+pow(y,2),2),atan2(x,y));
}

std::array<double,3> subtract_vectors(std::array<double,3> a,std::array<double,3> b)
{
    std::array<double,3> out;
    out[0]=a[0]-b[0];
    out[1]=a[1]-b[1];
    out[2]=a[2]-b[2];
    return out;
}

double dot_product(std::array<double,3> a,std::array<double,3> b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

typedef std::array<double,4> state_type;

struct push_back_state_and_time
{
    std::vector<state_type>& m_states;
    std::vector<double>& m_times;

    push_back_state_and_time(std::vector<state_type> &states,std::vector<double> &times)
    : m_states(states),m_times(times) {}

    void operator() (const state_type &x,double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

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
runge_kutta4<state_type> stepper;

class Ball
{
    public:
        double _x=0;
        double _y=0;
        double _z=ball_radius;
        double _r=ball_radius;
        bool _potted=false;
        bool _onrim=false;
        bool _inflight=false;

        //attributes for movement.
        double _rspin=0;
        double _xspin=0;
        double _yspin=0;
        double _vx=0;
        double _vy=0;
        double _vz=0;

        //grid position.
        int _gpos[4][2]={};
        int _order;

        //equation coefficients.
        std::array<double,3> _ax;
        std::array<double,3> _ay;
        std::array<double,3> _az;
        std::array<double,2> _avx;
        std::array<double,2> _avy;
        std::array<double,2> _avz;
        std::array<double,2> _awx;
        std::array<double,2> _awy;
        std::array<double,2> _awz;
        std::vector<std::array<double,3>> _rimpos;
        std::vector<std::array<double,3>> _rimvel;
        std::vector<double> _times;
        double _t;//till next phase change.

        //graphics.
        sf::CircleShape _shape;

        //class functions.
        Ball();
        void update_equation();
        std::array<double,3> get_position(double t);
        std::array<double,3> get_velocity(double t);
        std::array<double,3> get_spin(double t);
};

Ball::Ball()
{
    _shape.setRadius(_r*dfactor);
    _shape.setOrigin(dfactor*_r/2,dfactor*_r/2);
}

void Ball::update_equation()
{
    _ax[0]=_x;
    _ax[1]=_vx;
    _ax[2]=0.;
    _ay[0]=_y;
    _ay[1]=_vy;
    _ay[2]=0.;
    _az[0]=_z;
    _az[1]=_vz;
    _az[2]=0.;

    _avx[0]=_vx;
    _avx[1]=0.;
    _avy[0]=_vy;
    _avy[1]=0.;
    _avz[0]=_vz;
    _avz[1]=0.;

    _awx[0]=_xspin;
    _awx[1]=0.;
    _awy[0]=_yspin;
    _awy[1]=0.;
    _awz[0]=_rspin;
    _awz[1]=0.;

    _t=999.;

    double dist;

    _inflight=false;
    _onrim=false;

    double phi;
    double theta;
    double vphi;
    double vtheta;
    double wphi;
    double wtheta;

    double nspeed;
    double parspeed;
    double normal;

    double relx;
    double rely;

    std::array<double,3> temp;

    //check if on rim of pocket or falling within the pocket.

    for (int i=0;i<2;i++)
    {
        dist=nroot(pow(_x-mpockets[i][0],2)+pow(_y-mpockets[i][1],2),2);
        if (dist<mpocket_r-ball_radius)
        {
            //falling.
            _az[2]=-0.5*gravity;
            _avz[1]=-gravity;
            _inflight=true;
            break;
        }
        else if (dist>=mpocket_r-ball_radius && dist<mpocket_r+epsilon)
        {
            dist=nroot(pow(mpocket_r-dist,2)+pow(_z,2),2);
            if (dist<=ball_radius+epsilon)
            {
                //on the rim.
                //check if normal force is zero.
                phi=acos(_z/dist);
                nspeed=nroot(pow(_vx,2)+pow(_vy,2),2)*cos(atan2(mpockets[i][0]-_x,mpockets[i][1]-_y)-atan2(_vx,_vy));
                parspeed=nroot(pow(_vx,2)+pow(_vy,2),2)*cos(atan2(mpockets[i][0]-_x,mpockets[i][1]-_y)+0.5*pi-atan2(_vx,_vy));
                vtheta=parspeed/nroot(pow(_x-mpockets[i][0],2)+pow(_y-mpockets[i][1],2),2);
                vphi=nroot(pow(nspeed,2)+pow(_vz,2),2)/dist;
                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                //+ve parspeed means going clockwise round pocket.
                if (normal>0)
                {
                    //its on the rim m8.
                    _onrim=true;
                    wtheta=_xspin*sin(atan2(_x-mpockets[i][0],_y-mpockets[i][1]))+_yspin*cos(atan2(_x-mpockets[i][0],_y-mpockets[i][1]));
                    wphi=(_xspin*cos(atan2(_x-mpockets[i][0],_y-mpockets[i][1]))-_yspin*sin(atan2(_x-mpockets[i][0],_y-mpockets[i][1])))*cos(phi)-_rspin*sin(phi);
                    vtheta=(vtheta+(0.4*ball_radius*wphi)/(mpocket_r-ball_radius*sin(phi)))/(1.4);
                    vphi=(vphi-(0.4*ball_radius*wtheta)/(mpocket_r-ball_radius*sin(phi)))/(1.4);

                    double _angle=atan2(_vx,_vy);
                    relx=((_x-mpockets[i][0])-(_y-mpockets[i][1])*tan(_angle))/(cos(_angle)+sin(_angle)*tan(_angle));
                    rely=((_x-mpockets[i][0])*tan(_angle)+(_y-mpockets[i][1]))/(cos(_angle)+sin(_angle)*tan(_angle));
                    theta=pi-atan2(relx,rely);

                    //get the equation of motion set up.
                    _rimpos.clear();
                    _rimvel.clear();
                    _times.clear();
                    std::vector<state_type> output;
                    state_type xi={phi,vphi,theta,vtheta};

                    integrate_const(stepper,rimodem,xi,0.,0.5,0.001,push_back_state_and_time(output,_times));
                    size_t ind=0;

                    for (int j=0;j<500;j++)
                    {
                        //get upper limit on the time.
                        ind=j;
                        normal=gravity*cos(output[j][0])-ball_radius*pow(output[j][1],2)+(mpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][0])*pow(output[j][3],2);
                        //add on the xyz coords here.
                        relx=(mpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][2]);
                        rely=-(mpocket_r-ball_radius*sin(output[j][0]))*cos(output[j][2]);
                        temp={mpockets[i][0]+relx*cos(_angle)+rely*sin(_angle),mpockets[i][1]-relx*sin(_angle)+rely*cos(_angle),ball_radius*cos(output[j][0])};
                        _rimpos.push_back(temp);
                        std::tie(nspeed,theta)=add_vectors((mpocket_r-ball_radius*sin(output[j][0]))*output[j][3],atan2(mpockets[i][0]-temp[0],mpockets[i][1]-temp[1])+0.5*pi,ball_radius*cos(output[j][0])*output[j][1],atan2(mpockets[i][0]-temp[0],mpockets[i][1]-temp[1]));
                        temp={nspeed*sin(theta),nspeed*cos(theta),-ball_radius*sin(output[j][0])*output[j][1]};
                        _rimvel.push_back(temp);
                        if (normal<0)
                        {
                            break;
                        }
                    }
                    for (int j=0;j<500-ind;j++)
                    {
                        _times.pop_back();
                    }
                    _t=std::min(_times[_times.size()-1],_t);
                    break;
                }
                else
                {
                    _az[2]=-0.5*gravity;
                    _avz[1]=-gravity;
                    _inflight=true;
                    break;
                }
            }
            else
            {
                //falling.
                _az[2]=-0.5*gravity;
                _avz[1]=-gravity;
                _inflight=true;
                break;
            }
        }
    }
    for (int i=0;i<4;i++)
    {
        dist=nroot(pow(_x-cpockets[i][0],2)+pow(_y-cpockets[i][1],2),2);
        if (dist<cpocket_r-ball_radius)
        {
            //falling.
            _az[2]=-0.5*gravity;
            _avz[1]=-gravity;
            _inflight=true;
            break;
        }
        else if (dist>=cpocket_r-ball_radius && dist<cpocket_r+epsilon)
        {
            dist=nroot(pow(cpocket_r-dist,2)+pow(_z,2),2);
            if (dist<=ball_radius+epsilon)
            {
                //on the rim.
                //check if normal force is zero.
                phi=acos(_z/dist);
                nspeed=nroot(pow(_vx,2)+pow(_vy,2),2)*cos(atan2(cpockets[i][0]-_x,cpockets[i][1]-_y)-atan2(_vx,_vy));
                parspeed=nroot(pow(_vx,2)+pow(_vy,2),2)*cos(atan2(cpockets[i][0]-_x,cpockets[i][1]-_y)+0.5*pi-atan2(_vx,_vy));
                vtheta=parspeed/nroot(pow(_x-cpockets[i][0],2)+pow(_y-cpockets[i][1],2),2);
                vphi=nroot(pow(nspeed,2)+pow(_vz,2),2)/dist;
                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                //+ve parspeed means going clockwise round pocket.
                if (normal>0)
                {
                    //its on the rim m8.
                    _onrim=true;
                    wtheta=_xspin*sin(atan2(_x-cpockets[i][0],_y-cpockets[i][1]))+_yspin*cos(atan2(_x-cpockets[i][0],_y-cpockets[i][1]));
                    wphi=(_xspin*cos(atan2(_x-cpockets[i][0],_y-cpockets[i][1]))-_yspin*sin(atan2(_x-cpockets[i][0],_y-cpockets[i][1])))*cos(phi)-_rspin*sin(phi);
                    vtheta=(vtheta+(0.4*ball_radius*wphi)/(cpocket_r-ball_radius*sin(phi)))/(1.4);
                    vphi=(vphi-(0.4*ball_radius*wtheta)/(cpocket_r-ball_radius*sin(phi)))/(1.4);

                    double _angle=atan2(_vx,_vy);

                    relx=((_x-cpockets[i][0])-(_y-cpockets[i][1])*tan(_angle))/(cos(_angle)+sin(_angle)*tan(_angle));
                    rely=((_x-cpockets[i][0])*tan(_angle)+(_y-cpockets[i][1]))/(cos(_angle)+sin(_angle)*tan(_angle));
                    theta=pi-atan2(relx,rely);

                    //get the equation of motion set up.
                    _rimpos.clear();
                    _rimvel.clear();
                    _times.clear();
                    std::vector<state_type> output;
                    state_type xi={phi,vphi,theta,vtheta};

                    integrate_const(stepper,rimodec,xi,0.,0.5,0.001,push_back_state_and_time(output,_times));
                    size_t ind=0;

                    for (int j=0;j<500;j++)
                    {
                        //get upper limit on the time.
                        ind=j;
                        normal=gravity*cos(output[j][0])-ball_radius*pow(output[j][1],2)+(cpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][0])*pow(output[j][3],2);
                        //add on the xyz coords here.
                        relx=(cpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][2]);
                        rely=-(cpocket_r-ball_radius*sin(output[j][0]))*cos(output[j][2]);
                        temp={cpockets[i][0]+relx*cos(_angle)+rely*sin(_angle),cpockets[i][1]-relx*sin(_angle)+rely*cos(_angle),ball_radius*cos(output[j][0])};
                        _rimpos.push_back(temp);
                        std::tie(nspeed,theta)=add_vectors((cpocket_r-ball_radius*sin(output[j][0]))*output[j][3],atan2(cpockets[i][0]-temp[0],cpockets[i][1]-temp[1])+0.5*pi,ball_radius*cos(output[j][0])*output[j][1],atan2(cpockets[i][0]-temp[0],cpockets[i][1]-temp[1]));
                        temp={nspeed*sin(theta),nspeed*cos(theta),-ball_radius*sin(output[j][0])*output[j][1]};
                        _rimvel.push_back(temp);
                        if (normal<0)
                        {
                            break;
                        }
                    }
                    for (int j=0;j<500-ind;j++)
                    {
                        _times.pop_back();
                    }
                    _t=std::min(_times[_times.size()-1],_t);
                    break;
                }
                else
                {
                    _az[2]=-0.5*gravity;
                    _avz[1]=-gravity;
                    _inflight=true;
                    break;
                }
            }
            else
            {
                //falling.
                _az[2]=-0.5*gravity;
                _avz[1]=-gravity;
                _inflight=true;
                break;
            }
        }
    }

    if (!_onrim)
    {
        //check if in ballistic motion.
        //pockets first.
        if (_z>ball_radius || _vz>0)
        {
            _inflight=true;
            _az[2]=-0.5*gravity;
            _avz[1]=-gravity;
        }

        //on bed of table and not moving up or down.
        if (!_inflight)
        {
            _awz[1]=-sgn(_rspin)*mus*gravity;
            double m=nroot(pow(_vx-_r*_xspin,2)+pow(_vy-_r*_yspin,2),2);
            if ((fabs(_xspin*_r-_vx)<epsilon && fabs(_yspin*_r-_vy)<epsilon) || fabs(m)<pow(10,-15))
            {
                //both rolling.
                m=nroot(pow(_vx,2)+pow(_vy,2),2);
            }

            if (fabs(_vx)<epsilon && fabs(_xspin)<epsilon)
            {
                //stopped in x.
                _ax[2]=0;
                _avx[1]=0;
                _awx[1]=0;
            }
            else if (fabs(_xspin*_r-_vx)<epsilon)
            {
                //rolling in x.
                _vx=_xspin*_r;
                _ax[2]=-0.5*mus*gravity*_vx/m;
                _t=std::min(m/(mus*gravity),_t);
                _avx[1]=-mus*gravity*_vx/m;
                _awx[1]=-mus*gravity*_vx/(_r*m);
            }
            else
            {
                //sliding in x.
                _ax[2]=0.5*muk*gravity*(_xspin*_r-_vx)/m;
                _t=std::min((2./7.)*m/(muk*gravity),_t);
                _avx[1]=(_xspin*_r-_vx)*muk*gravity/m;
                _awx[1]=-(_xspin*_r-_vx)*2.5*muk*gravity/(_r*m);
            }

            if (fabs(_vy)<epsilon && fabs(_yspin)<epsilon)
            {
                //stopped in y.
                _ay[2]=0;
                _avy[1]=0;
                _awy[1]=0;
            }
            else if (fabs(_yspin*_r-_vy)<epsilon)
            {
                //rolling in y.
                _vy=_yspin*_r;
                _ay[2]=-0.5*mus*gravity*_vy/m;
                _t=std::min(m/(mus*gravity),_t);
                _avy[1]=-mus*gravity*_vy/m;
                _awy[1]=-mus*gravity*_vy/(_r*m);
            }
            else
            {
                //sliding in y.
                _ay[2]=0.5*muk*gravity*(_yspin*_r-_vy)/m;
                _t=std::min((2./7.)*m/(muk*gravity),_t);
                _avy[1]=(_yspin*_r-_vy)*muk*gravity/m;
                _awy[1]=-(_yspin*_r-_vy)*2.5*muk*gravity/(_r*m);
            }
        }
    }

    if (_inflight)
    {
        std::array<double,2> roots=qsolve_quadratic(_az[2],_az[1],_az[0]+ball_radius);
        for (int i=0;i<2;i++)
        {
            if (roots[i]==roots[i] && roots[i]>=epsilon && roots[i]<_t)
            {
                phi=roots[i];
                theta=_az[2]*phi*phi+_az[1]*phi+_az[0]+ball_radius;
                if (theta<0)
                {
                    _t=phi;
                }
                else
                {
                    int c=0;
                    while (true)
                    {
                        c+=1;
                        phi=roots[i]-c*epsilon;
                        theta=_az[2]*phi*phi+_az[1]*phi+_az[0]+ball_radius;
                        if (theta<0 && phi>=epsilon)
                        {
                            _t=phi;
                            break;
                        }
                        phi=roots[i]+c*epsilon;
                        theta=_az[2]*phi*phi+_az[1]*phi+_az[0]+ball_radius;
                        if (theta<0 && phi<_t)
                        {
                            _t=phi;
                            break;
                        }
                        if (c>1000)
                        {
                            break;
                        }
                    }
                }
            }
        }
    }

    if (_order==1 && _onrim)
    {
        std::cout << "On the rim!" << std::endl;
        for (int i=0;i<1;i++)
        {
            std::cout << _rimpos[i][0] << " : " << _rimpos[i][1] << " : " << _rimpos[i][2] << " : " << _times[i] << std::endl;
            std::cout << _rimvel[i][0] << " : " << _rimvel[i][1] << " : " << _rimvel[i][2] << std::endl;
        }
    }
    if (_order==1 && _inflight)
    {
        std::cout << "In flight!" << std::endl;
    }
//    if (_order==1 && !_inflight && !_onrim && !_potted)
//    {
//        std::cout << "On bed!" << std::endl;
//    }
}

std::array<double,3> Ball::get_position(double t)
{
    std::array<double,3> pos;
    pos[0]=_ax[0]+_ax[1]*t+_ax[2]*pow(t,2);
    pos[1]=_ay[0]+_ay[1]*t+_ay[2]*pow(t,2);
    pos[2]=_az[0]+_az[1]*t+_az[2]*pow(t,2);
    return pos;
}

std::array<double,3> Ball::get_velocity(double t)
{
    std::array<double,3> vel;
    vel[0]=_avx[0]+_avx[1]*t;
    vel[1]=_avy[0]+_avy[1]*t;
    vel[2]=_avz[0]+_avz[1]*t;
    return vel;
}

std::array<double,3> Ball::get_spin(double t)
{
    std::array<double,3> spin;
    spin[0]=_awx[0]+_awx[1]*t;
    spin[1]=_awy[0]+_awy[1]*t;
    spin[2]=_awz[0]+_awz[1]*t;
    return spin;
}

class Cushion
{
    public:
        double _x;
        double _y;
        double _angle;
        double _length;
        int _p1;
        int _p2;

        sf::ConvexShape _shape;
        sf::ConvexShape _railshape;
        sf::ConvexShape _pocketshape;
        sf::ConvexShape _p1shape;
        sf::ConvexShape _p2shape;

        Cushion();
        Cushion(double x1,double y1,double angle1,int p1,int p2);
        std::tuple<double,double> distance(double x,double y,double height);
};

Cushion::Cushion(){}

Cushion::Cushion(double x1,double y1,double angle1,int p1,int p2)
{
    _x=x1;
    _y=y1;
    _angle=angle1;
    _p1=p1;
    _p2=p2;

    _railshape.setFillColor(sf::Color(railcolour[0],railcolour[1],railcolour[2]));
    _railshape.setPointCount(6);

    _p1shape.setFillColor(sf::Color(leathercolour[0],leathercolour[1],leathercolour[2]));
    _p2shape.setFillColor(sf::Color(leathercolour[0],leathercolour[1],leathercolour[2]));

    _pocketshape.setFillColor(sf::Color(0,0,0));

    _shape.setFillColor(sf::Color(cushioncolour[0],cushioncolour[1],cushioncolour[2]));

    double points[6][2];
    double k1=8.595-9*sin(pi/4);
    double delta=0.5;
    double x;
    double y;
    double angle;

    if (_p1!=0 || _p2!=0)
    {
        _length=cush_thickness+table_width;

        if (_p1==0)
        {
            points[0][0]=k1;
            points[0][1]=-2.445;
            points[1][0]=k1;
            points[1][1]=0;
            points[2][0]=k1+delta;
            points[2][1]=0;
            points[3][0]=_length-1.625-0.438;
            points[3][1]=0;
            points[4][0]=_length-1.625;
            points[4][1]=-0.438;
            points[5][0]=_length-1.625;
            points[5][1]=-2.445;
        }
        else
        {
            points[0][0]=1.625;
            points[0][1]=-2.445;
            points[1][0]=1.625;
            points[1][1]=-0.438;
            points[2][0]=1.625+0.438;
            points[2][1]=0;
            points[3][0]=_length-(k1+delta);
            points[3][1]=0;
            points[4][0]=_length-k1;
            points[4][1]=0;
            points[5][0]=_length-k1;
            points[5][1]=-2.445;
        }
    }
    else
    {
        _length=2*cush_thickness+table_width;

        points[0][0]=k1;
        points[0][1]=-2.445;
        points[1][0]=k1;
        points[1][1]=0;
        points[2][0]=k1+delta;
        points[2][1]=0;
        points[3][0]=_length-(k1+delta);
        points[3][1]=0;
        points[4][0]=_length-k1;
        points[4][1]=0;
        points[5][0]=_length-k1;
        points[5][1]=-2.445;
    }

    for (int i=0;i<6;i++)
    {
        x=points[i][0];
        y=points[i][1];
        _railshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
    }

    //pocket hole.
    if (_p1==0)
    {
        _pocketshape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=0.5*pi*i/50;
            x=cpocket_r*sin(angle);
            y=cpocket_r*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=0.75*pi+pi*(i-51)/50;
            x=0.5*k1+0.5*k1*sqrt(2.0)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.0)*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }
    else
    {
        _pocketshape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=-0.5*pi+pi*i/50;
            x=mpocket_r*sin(angle);
            y=-0.156+mpocket_r*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi-asin(1.625/pround)+2*asin(1.625/pround)*(i-51)/50;
            x=pround*sin(angle);
            y=1.875+pround*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }

    //pocket leathers.
    if (_p1==0)
    {
        _p1shape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=pi*1.25-0.5*pi*i/50;
            x=0.5*k1+0.5*k1*sqrt(2.0)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.0)*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi+0.25*pi*(i-51)/50;
            x=k1+(k1+2.445)*sin(angle);
            y=k1+(k1+2.445)*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }
    else
    {
        _p1shape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=pi-asin(1.625/pround)*i/50;
            x=pround*sin(angle);
            y=1.875+pround*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi-asin(1.625/pround)+asin(1.625/pround)*(i-51)/50;
            x=pround*sin(angle);
            y=1.875-(4.32-pround)+pround*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }
    if (_p2==0)
    {
        _p2shape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=1.25*pi-0.5*pi*i/50;
            x=_length-0.5*k1+0.5*k1*sqrt(2.0)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.0)*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=0.75*pi+0.25*pi*(i-51)/50;
            x=_length-k1+(k1+2.445)*sin(angle);
            y=k1+(k1+2.445)*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }
    else
    {
        _p2shape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=pi+asin(1.625/pround)-asin(1.625/pround)*i/50;
            x=_length+pround*sin(angle);
            y=1.875+pround*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi+asin(1.625/pround)*(i-51)/50;
            x=_length+pround*sin(angle);
            y=1.875-(4.32-pround)+pround*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }

    //get shapes for pockets.
    int pos=0;

    _shape.setPointCount(202);

    if (p1==0)
    {
        for (int i=50;i<51;i++)
        {
            angle=pi*5/4-i*(pi/2)/50;
            x=0.5*k1+0.5*k1*sqrt(2.0)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.0)*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=pi*7/4+(pi/4-atan(6/11))*i/50;
            x=6.15+4.5*sin(angle);
            y=-2.445+4.5*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=atan(6/11)*(i/50-1);
            x=5.43+3*sin(angle);
            y=-1.125+3*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
    }
    else
    {
//        for (int i=0;i<50;i++)
//        {
//            angle=pi*2-(pi/4)*i/50;
//            x=1.625*sin(angle);
//            y=-0.438+1.625*cos(angle);
//            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
//            pos+=1;
//        }
        for (int i=0;i<50;i++)
        {
            angle=pi*3/2+atan(0.812/0.625)*i/50;
            x=3.719+2.094*sin(angle);
            y=-0.438+2.094*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<51;i++)
        {
            angle=pi*3/2+atan(0.812/0.625)+atan(0.625/0.812)*i/50;
            x=4.344+3.125*sin(angle);
            y=-1.25+3.125*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
    }
    if (p2==0)
    {
        for (int i=0;i<50;i++)
        {
            angle=atan(6/11)*i/50;
            x=_length-5.43+3*sin(angle);
            y=-1.125+3*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=atan(6/11)+(pi/4-atan(6/11))*i/50;
            x=_length-6.15+4.5*sin(angle);
            y=-2.445+4.5*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<1;i++)
        {
            angle=pi*5/4-i*(pi/2)/50;
            x=_length-0.5*k1+0.5*k1*sqrt(2.0)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.0)*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
    }
    else
    {
        for (int i=0;i<50;i++)
        {
            angle=atan(0.625/0.812)*i/50;
            x=_length-4.344+3.125*sin(angle);
            y=-1.25+3.125*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<51;i++)
        {
            angle=atan(0.625/0.812)+(pi/2-atan(0.625/0.812))*i/50;
            x=_length-3.719+2.094*sin(angle);
            y=-0.438+2.094*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
//        for (int i=0;i<50;i++)
//        {
//            angle=pi*2-(pi/4)*i/50;
//            x=_length-(1.625*sin(angle));
//            y=-0.438+1.625*cos(angle);
//            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
//            pos+=1;
//        }
    }
}

std::tuple<double,double> Cushion::distance(double x,double y,double height)
{
    //find closest distance from a point to the cushion and return angle of normal vector.
    double relx=((x-_x)-(y-_y)*tan(_angle))/(cos(_angle)+sin(_angle)*tan(_angle));
    double rely=((x-_x)*tan(_angle)+(y-_y))/(cos(_angle)+sin(_angle)*tan(_angle));

    double min_dist=999;
    double normal_angle;

    double a1;
    double a2;

    double dist;
    double angle;

    //check pocket closest to the rel origin.
    //corner pocket.
    if (_p1==0)
    {
        if (height<0)
        {
            angle=atan2(relx,rely);
            if (angle>=0.25*pi && angle<cpocket_angle)
            {
                dist=nroot(pow(relx,2)+pow(rely,2),2);
                if (dist<cpocket_r)
                {
                    dist=cpocket_r-dist;
                    if (dist<min_dist)
                    {
                        min_dist=dist;
                        normal_angle=angle+pi;
                    }
                }
            }
        }
        //check in pocket first.
        angle=atan2(k1*0.5-relx,k1*0.5-rely);
        if (angle>-pi/4. && angle<=pi/4.)
        {
            dist=fabs(nroot(pow(k1*0.5-relx,2)+pow(k1*0.5-rely,2),2)-0.5*k1*nroot(2.,2));
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check straight bit.
        a1=0.5*(relx+rely+k1);
        a2=a1-k1;
        if (a1>=k1 && a1<k2)
        {
            dist=nroot(pow(a1-relx,2)+pow(a2-rely,2),2);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=pi*7./4.;
            }
        }
        //check curvy 1.
        angle=atan2(relx-6.15,rely+2.445);
        if (angle>=-pi/4. && angle<-atan(6./11.))
        {
            dist=fabs(nroot(pow(relx-6.15,2)+pow(rely+2.445,2),2)-4.5);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check curvy 2.
        angle=atan2(relx-5.43,rely+1.125);
        if (angle>=-atan(6./11.) && angle<0)
        {
            dist=fabs(nroot(pow(relx-5.43,2)+pow(rely+1.125,2),2)-3.);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //middle pocket next.
        if (_p2!=0)
        {
            //check straight bit.
            if (relx>=5.43 && relx<_length-4.344)
            {
                dist=fabs(rely-cush_thickness);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=0;
                }
            }
            //check curvy.
            angle=atan2(relx-_length+4.344,rely+1.25);
            if (angle>=0 && angle<atan(0.625/0.812))
            {
                dist=fabs(nroot(pow(relx-_length+4.344,2)+pow(rely+1.25,2),2)-3.125);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            //check curvy2.
            angle=atan2(relx-_length+3.719,rely+0.438);
            if (angle>=atan(0.625/0.812) && angle<pi/2.)
            {
                dist=fabs(nroot(pow(relx-_length+3.719,2)+pow(rely+0.438,2),2)-2.094);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            //check straight bit.
            if (rely<=-0.438 && rely>(1.875-pround*cos(asin(1.625/pround))))
            {
                dist=fabs(_length-1.625-relx);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=0.5*pi;
                }
            }
            //check in pocket.
            angle=atan2(-relx+_length,-rely+1.875);
            if (angle<=asin(1.625/pround) && angle>0)
            {
                dist=fabs(nroot(pow(relx-_length,2)+pow(-rely+1.875,2),2)-pround);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
        }
        //corner pocket next.
        else
        {
            //check straight bit.
            if (relx>=5.43 && relx<_length-5.43)
            {
                dist=fabs(rely-cush_thickness);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=0;
                }
            }
            //check curvy.
            angle=atan2(relx-_length+5.43,rely+1.125);
            if (angle>=0 && angle<atan(6./11.))
            {
                dist=fabs(nroot(pow(relx-_length+5.43,2)+pow(rely+1.125,2),2)-3.);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            //check curvy 2.
            angle=atan2(relx-_length+6.15,rely+2.445);
            if (angle>=atan(6./11.) && angle<pi/4.)
            {
                dist=fabs(nroot(pow(relx-_length+6.15,2)+pow(rely+2.445,2),2)-4.5);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            //check straight bit.
            a1=0.5*(_length-k1+relx-rely);
            a2=a1+rely-relx;
            if (a1>=_length-k2 && a1<_length-k1)
            {
                dist=nroot(pow(relx-a1,2)+pow(rely-a2,2),2);
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=pi/4.;
                }
            }
            //check in pocket.
            angle=atan2(_length-k1*0.5-relx,k1*0.5-rely);
            if (angle<=pi/4. && angle>-pi/4.)
            {
                dist=fabs(nroot(pow(k1*0.5-rely,2)+pow(_length-k1*0.5-relx,2),2)-0.5*k1*nroot(2.,2));
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            if (height<0)
            {
                angle=atan2(relx-_length,rely);
                if (angle>=-cpocket_angle && angle<-0.25*pi)
                {
                    dist=nroot(pow(relx-_length,2)+pow(rely,2),2);
                    if (dist<cpocket_r)
                    {
                        dist=cpocket_r-dist;
                        if (dist<min_dist)
                        {
                            min_dist=dist;
                            normal_angle=angle+pi;
                        }
                    }
                }
            }
        }
    }
    //middle pocket.
    else
    {
        //check in pocket first.
        angle=atan2(-relx,-rely+1.875);
        if (angle<=0 && angle>-asin(1.625/pround))
        {
            dist=fabs(nroot(pow(-rely+1.875,2)+pow(relx,2),2)-pround);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check straight bit.
        if (rely<=-0.438 && rely>(1.875-pround*cos(asin(1.625/pround))))
        {
            dist=fabs(relx-1.625);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=1.5*pi;
            }
        }
        //check curvy.
        angle=atan2(relx-3.719,rely+0.438);
        if (angle>=-pi/2 && angle<-atan(0.625/0.812))
        {
            dist=fabs(nroot(pow(rely+0.438,2)+pow(relx-3.719,2),2)-2.094);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check curvy 2.
        angle=atan2(relx-4.344,rely+1.25);
        if (angle>=-atan(0.625/0.812) && angle<0)
        {
            dist=fabs(nroot(pow(relx-4.344,2)+pow(rely+1.25,2),2)-3.125);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check straight bit.
        if (relx>=4.344 && relx<_length-5.43)
        {
            dist=fabs(rely-cush_thickness);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=0;
            }
        }
        //check curvy.
        angle=atan2(relx-_length+5.43,rely+1.125);
        if (angle>=0 && angle<atan(6./11.))
        {
            dist=fabs(nroot(pow(relx-_length+5.43,2)+pow(rely+1.125,2),2)-3.);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check curvy 2.
        angle=atan2(relx-_length+6.15,rely+2.445);
        if (angle>=atan(6./11.) && angle<pi/4.)
        {
            dist=fabs(nroot(pow(relx-_length+6.15,2)+pow(rely+2.445,2),2)-4.5);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check straight bit.
        a1=0.5*(_length-k1+relx-rely);
        a2=a1+rely-relx;
        if (a1>=_length-k2 && a1<_length-k1)
        {
            dist=nroot(pow(relx-a1,2)+pow(rely-a2,2),2);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=pi/4.;
            }
        }
        //check in pocket.
        angle=atan2(_length-k1*0.5-relx,k1*0.5-rely);
        if (angle<=pi/4. && angle>-pi/4.)
        {
            dist=fabs(nroot(pow(k1*0.5-rely,2)+pow(_length-k1*0.5-relx,2),2.)-0.5*k1*nroot(2.,2));
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
    }

    //rotate the normal angle.
    normal_angle=normal_angle+_angle;

    return std::make_tuple(min_dist,normal_angle);
}

class Cue
{
    public:
        double _angle=-0.5*pi;

        sf::Texture _texture;
        sf::Sprite _sprite;

        Cue();
};

Cue::Cue()
{
    if (!_texture.loadFromFile("cue.png"))
    {
        std::cout << "Error loading cue texture!" << std::endl;
    }
    _texture.setSmooth(true);
    _sprite.setTexture(_texture);
    _sprite.scale(58.*dfactor/1369.,58.*dfactor/1369.);
    _sprite.setOrigin(0,26.*58.*dfactor/1369.);
}

Eigen::Matrix<double,46,1> collisions(Ball b[22],Cushion cush[6])
{
    //broad phase collision check.
    Eigen::MatrixXd dv=Eigen::MatrixXd::Zero(46,1);

//    for (int i=0;i<46;i++)
//    {
//        dv(i,0)=0;
//    }

    int val;
    int thing;
    double dx;
    double dy;
    double dist;
    double relspeed;

    std::set<std::tuple<int,int>> collided;
    std::vector<int> already;

    Eigen::MatrixXd Z(46,0);
    int col=0;

    for (int i=0;i<22;i++)
    {
        if (b[i]._potted==false)
        {
            if (b[i]._vx!=0 || b[i]._vy!=0 || b[i]._vz!=0)
            {
                already.clear();
                for (int j=0;j<4;j++)
                {
                    val=grid_index[b[i]._gpos[j][0]][b[i]._gpos[j][1]];
                    if (val>1)
                    {
                        for (int k=0;k<val;k++)
                        {
                            //add to potentially touching list.
                            thing=grid[b[i]._gpos[j][0]][b[i]._gpos[j][1]][k];
                            if (thing!=b[i]._order && b[thing-1]._potted==false)
                            {
                                if (std::find(already.begin(),already.end(),thing)==already.end())
                                {
                                    if (i>(thing-1) && collided.find(std::make_tuple(thing-1,i))!=collided.end())
                                    {
                                        continue;
                                    }
                                    else if (i<(thing-1) && collided.find(std::make_tuple(i,thing-1))!=collided.end())
                                    {
                                        continue;
                                    }
                                    already.push_back(thing);
                                    //get the distance using pythag.
                                    dx=b[thing-1]._x-b[i]._x;
                                    dy=b[thing-1]._y-b[i]._y;
                                    dist=nroot(pow(dx,2)+pow(dy,2),2.);
                                    relspeed=nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2.)*cos(atan2(dx,dy)-atan2(b[i]._vx,b[i]._vy))-nroot(pow(b[thing-1]._vx,2)+pow(b[thing-1]._vy,2),2.)*cos(atan2(dx,dy)-atan2(b[thing-1]._vx,b[thing-1]._vy));
                                    if (dist-2.*ball_radius<epsilon && relspeed>0)
                                    {
                                        //actual collision.
                                        if (i>(thing-1))
                                        {
                                            collided.insert(std::make_tuple(thing-1,i));
                                        }
                                        else
                                        {
                                            collided.insert(std::make_tuple(i,thing-1));
                                        }
                                        Z.conservativeResizeLike(Eigen::MatrixXd::Zero(46,col+1));
                                        Z(2*i,col)=-dx/dist;
                                        Z(2*i+1,col)=-dy/dist;
                                        Z(2*(thing-1),col)=dx/dist;
                                        Z(2*(thing-1)+1,col)=dy/dist;
                                        std::cout << i << ":" << thing-1 << std::endl;
                                        col+=1;
                                    }
                                }
                            }
                        }
                    }
                }
                //check cushion collisions.
                if (fabs(blue_x-b[i]._x)>table_width-ball_radius-100*epsilon || fabs(blue_y-b[i]._y)>0.5*table_width-ball_radius-100*epsilon)
                {
                    double min_dist=999;
                    double normal_angle;
                    double angle;
                    for (int j=0;j<6;j++)
                    {
                        std::tie(dist,angle)=cush[j].distance(b[i]._x,b[i]._y,b[i]._z);
                        if (dist<min_dist)
                        {
                            min_dist=dist;
                            normal_angle=angle;
                        }
                    }
                    relspeed=-nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*cos(normal_angle-atan2(b[i]._vx,b[i]._vy));
                    if (min_dist-ball_radius<100*epsilon && relspeed>0)
                    {
                        //collision with cushion.
                        Z.conservativeResizeLike(Eigen::MatrixXd::Zero(46,col+1));
                        Z(2*i,col)=sin(normal_angle);
                        Z(2*i+1,col)=cos(normal_angle);
                        Z(44,col)=-sin(normal_angle);
                        Z(45,col)=-cos(normal_angle);
                        col+=1;
                        //adjust the vertical spin off cushion.
                        double relspin=b[i]._xspin*sin(normal_angle+pi)+b[i]._yspin*cos(normal_angle+pi);
                        double dw;
                        if (relspin>0)
                        {
                            //topspin.
                            dw=-5.*relspeed*(muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2);
                        }
                        if (relspin<0)
                        {
                            //backspin.
                            dw=-5.*relspeed*(-muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2);
                        }
                        if (relspin==0)
                        {
                            dw=-5.*relspeed*cushion_diff/pow(ball_radius,2);
                        }
                        dw=fmax(dw,-relspin-relspeed/ball_radius);
                        b[i]._xspin+=dw*sin(normal_angle+pi);
                        b[i]._yspin+=dw*cos(normal_angle+pi);
                        //adjust sidespin off cushion.
                        double parspeed=nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*sin(atan2(b[i]._vx,b[i]._vy)-normal_angle);
                        double dvpar;

                        if ((parspeed>=0 && b[i]._rspin*ball_radius>parspeed) || (parspeed<0 && b[i]._rspin*ball_radius>parspeed))
                        {
                            dw=-5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius<parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        else if ((parspeed>=0 && b[i]._rspin*ball_radius<parspeed) || (parspeed<0 && b[i]._rspin*ball_radius<parspeed))
                        {
                            dw=5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius>parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        b[i]._rspin+=dw;
                        double _speed;
                        double _angle;
                        std::tie(_speed,_angle)=add_vectors(nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2),atan2(b[i]._vx,b[i]._vy),dvpar,normal_angle+0.5*pi);
                        b[i]._vx=_speed*sin(_angle);
                        b[i]._vy=_speed*cos(_angle);
                    }
                }
            }
        }
    }

    //check if any total collisions.
    if (col==1)
    {
        //single collision.
        //check if to cushion.
        if (Z(44,0)!=0 || Z(45,0)!=0)
        {
            double v;
            double an;
            for (int i=0;i<22;i++)
            {
                if (Z(2*i,0)!=0 || Z(2*i+1,0)!=0)
                {
                    //collision with this ball.
                    std::tie(v,an)=add_vectors(-nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*cos(atan2(b[i]._vx,b[i]._vy)-atan2(Z(2*i,0),Z(2*i+1,0))),atan2(Z(2*i,0),Z(2*i+1,0)),nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*cos(atan2(b[i]._vx,b[i]._vy)-(atan2(Z(2*i,0),Z(2*i+1,0))+0.5*pi)),atan2(Z(2*i,0),Z(2*i+1,0))+0.5*pi);
                    dv(2*i)=v*sin(an)-b[i]._vx;
                    dv(2*i+1)=v*cos(an)-b[i]._vy;
                    break;
                }
            }
        }
        //ball-ball collision.
        else
        {
            for (int i=0;i<22;i++)
            {
                if (Z(2*i,0)!=0 || Z(2*i+1,0)!=0)
                {
                    for (int j=i+1;j<22;j++)
                    {
                        if (Z(2*j,0)!=0 || Z(2*j+1,0)!=0)
                        {
                            double na=atan2(Z(2*i,0),Z(2*i+1,0));
                            double pari=nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*cos(atan2(b[i]._vx,b[i]._vy)-na);
                            double peri=nroot(pow(b[i]._vx,2)+pow(b[i]._vy,2),2)*cos(atan2(b[i]._vx,b[i]._vy)-(na+0.5*pi));
                            double parj=nroot(pow(b[j]._vx,2)+pow(b[j]._vy,2),2)*cos(atan2(b[j]._vx,b[j]._vy)-na);
                            double perj=nroot(pow(b[j]._vx,2)+pow(b[j]._vy,2),2)*cos(atan2(b[j]._vx,b[j]._vy)-(na+0.5*pi));
                            double speed;
                            double theta;
                            std::tie(speed,theta)=add_vectors(parj,na,peri,na+0.5*pi);
                            dv(2*i)=speed*sin(theta)-b[i]._vx;
                            dv(2*i+1)=speed*cos(theta)-b[i]._vy;
                            std::tie(speed,theta)=add_vectors(pari,na,perj,na+0.5*pi);
                            dv(2*j)=speed*sin(theta)-b[j]._vx;
                            dv(2*j+1)=speed*cos(theta)-b[j]._vy;
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    //simultaneous.
    if (col>1)
    {
        Eigen::MatrixXd v(46,1);
        Eigen::MatrixXd J(col,1);

        for (int i=0;i<22;i++)
        {
            //v(2*i,0)=b[i]._speed*sin(b[i]._angle);
            //v(2*i+1,0)=b[i]._speed*cos(b[i]._angle);
            v(2*i,0)=b[i]._vx;
            v(2*i+1,0)=b[i]._vy;
        }
        v(44,0)=0;
        v(45,0)=0;

        //J=(Z.transpose()*M_*Z).colPivHouseholderQr().solve(-(1.+eball)*Z.transpose()*v);
        J=(Z.transpose()*M_*Z).householderQr().solve(-(1.+eball)*Z.transpose()*v);
        dv=M_*Z*J;

//        std::cout << "Z" << std::endl;
//        std::cout << Z << std::endl;
//        std::cout << "J" << std::endl;
//        std::cout << J << std::endl;
//        std::cout << "dv" << std::endl;
//        std::cout << dv << std::endl;
    }
    return dv;
}

std::vector<std::array<double,66>> simulate(Ball balls[22],Cushion cush[6])
{
    //assumes that the input is not a static scenario (where all balls are still).
    std::vector<std::array<double,66>> pos;
    std::array<double,66> temp;
    std::array<double,3> out;
    std::array<double,3> out2;
    std::array<double,3> out3;
    std::array<double,3> out4;
    std::array<double,4> quartic;
    std::array<double,2> quadratic;
    std::array<std::array<double,5>,6> ctemp;
    std::array<double,5> x2;
    std::array<double,5> y2;
    std::array<double,5> z2;
    std::array<double,9> toctic;
    std::array<double,8> monoctic;
    std::array<double,8> octic;
    Eigen::Matrix<double,46,1> dv;
    Eigen::Matrix<double,46,1> zero;

    for (int i=0;i<46;i++)
    {
        dv(i,0)=0;
        zero(i,0)=0;
    }

    //set up the equations for each ball.

    double t; //time until next event.
    double totaltime=0; //total time elapsed.
    double start;
    int c;

    //for solving quartics.
    double a1;
    double b1;
    double c1;
    double a2;
    double b2;
    double c2;
    double a3;
    double b3;
    double c3;
    double t0;
    double t1;
    double t2;
    double t3;
    double t4;
    double x;
    double y;
    double phi;
    double nspeed;
    double parspeed;
    double vtheta;
    double vphi;
    double normal;

    double xmin;
    double xmax;
    double xmin2;
    double xmax2;
    double ymin;
    double ymax;
    double ymin2;
    double ymax2;

    int xpos;
    int ypos;
    int gxmin;
    int gxmax;
    int gymin;
    int gymax;

    bool still=false;
    bool collision=false;

    while (!still)
    {
        std::cout << "Time: " << totaltime << std::endl;

        collision=false;
        t=999.;
        still=true;
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted==false)
            {
                balls[i].update_equation();
                if (balls[i]._vx!=0 || balls[i]._vy!=0 || balls[i]._vz!=0 || balls[i]._xspin!=0 || balls[i]._yspin!=0)
                {
                    still=false;
                }
                if (balls[i]._t<t)
                {
                    t=balls[i]._t;
                }
            }
        }
        if (still)
        {
            break;
        }
        //check for collisions (cushion/balls) or if ball into pocket.
        //check moving balls only.
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted || balls[i]._onrim)
            {
                continue;
            }
            if (balls[i]._vx==0 && balls[i]._vy==0 && balls[i]._vz==0 && balls[i]._xspin==0 && balls[i]._yspin==0)
            {
                continue;
            }

            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4*balls[i]._ax[2]));
            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4*balls[i]._ax[2]));
            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4*balls[i]._ay[2]));
            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4*balls[i]._ay[2]));
            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);

            //check for ball-ball collisions.
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted || balls[j]._onrim)
                {
                    continue;
                }
                if (i==j)
                {
                    continue;
                }

                //initial check for i.
                xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4*balls[j]._ax[2]));
                xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4*balls[j]._ax[2]));
                xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4*balls[j]._ay[2]));
                ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4*balls[j]._ay[2]));
                ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);

                if (!((xmax>xmin2-2*ball_radius || xmax2>xmin-2*ball_radius) && (ymax>ymin2-2*ball_radius || ymax2>ymin-2*ball_radius)))
                {
                    continue;
                }
                //check if their trajectories are ever 2r apart.
                a1=balls[i]._ax[2]-balls[j]._ax[2];
                b1=balls[i]._ax[1]-balls[j]._ax[1];
                c1=balls[i]._ax[0]-balls[j]._ax[0];
                a2=balls[i]._ay[2]-balls[j]._ay[2];
                b2=balls[i]._ay[1]-balls[j]._ay[1];
                c2=balls[i]._ay[0]-balls[j]._ay[0];
                a3=balls[i]._az[2]-balls[j]._az[2];
                b3=balls[i]._az[1]-balls[j]._az[1];
                c3=balls[i]._az[0]-balls[j]._az[0];

                t0=pow(c1,2)+pow(c2,2)+pow(c3,2)-4*pow(ball_radius,2);
                t1=2.*(b1*c1+b2*c2+b3*c3);
                t2=2.*(a1*c1+a2*c2+a3*c3)+pow(b1,2)+pow(b2,2)+pow(b3,2);
                t3=2.*(a1*b1+a2*b2+a3*b3);
                t4=pow(a1,2)+pow(a2,2)+pow(a3,2);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=epsilon && quartic[k]<t)
                    {
                        //verify the time of actual collision to be certain.
                        t0=quartic[k];
                        out=balls[i].get_position(t0);
                        out2=balls[j].get_position(t0);
                        out3=balls[i].get_velocity(t0);
                        out4=balls[j].get_velocity(t0);
                        //a1=nroot(pow(out3[0],2)+pow(out3[1],2),2)*cos(atan2(out3[0],out3[1])-atan2(out2[0]-out[0],out2[1]-out[1]))+nroot(pow(out4[0],2)+pow(out4[1],2),2)*cos(atan2(out4[0],out4[1])-atan2(out[0]-out2[0],out[1]-out2[1]));
                        a1=dot_product(out3,subtract_vectors(out2,out))+dot_product(out4,subtract_vectors(out,out2));
                        if (nroot(pow(out[0]-out2[0],2)+pow(out[1]-out2[1],2)+pow(out[2]-out2[2],2),2)-2*ball_radius<epsilon && a1>0)
                        {
                            //collision at the specified time.
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            //adjust the time minutely to ensure a collision.
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[j].get_position(t0-c*epsilon);
                                out3=balls[i].get_velocity(t0-c*epsilon);
                                out4=balls[j].get_velocity(t0-c*epsilon);
                                //a1=nroot(pow(out3[0],2)+pow(out3[1],2),2)*cos(atan2(out3[0],out3[1])-atan2(out2[0]-out[0],out2[1]-out[1]))+nroot(pow(out4[0],2)+pow(out4[1],2),2)*cos(atan2(out4[0],out4[1])-atan2(out[0]-out2[0],out[1]-out2[1]));
                                a1=dot_product(out3,subtract_vectors(out2,out))+dot_product(out4,subtract_vectors(out,out2));
                                if (nroot(pow(out[0]-out2[0],2)+pow(out[1]-out2[1],2)+pow(out[2]-out2[2],2),2.)-2*ball_radius<epsilon && a1>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[j].get_position(t0+c*epsilon);
                                out3=balls[i].get_velocity(t0+c*epsilon);
                                out4=balls[j].get_velocity(t0+c*epsilon);
                                //a1=out3[0]*(out2[0]-out[0])+out3[1]*(out2[1]-out[1])+out3[2]*(out2[2]-out[2])+out4[0]*(out[0]-out2[0])+out4[1]*(out[1]-out2[1])+out4[2]*(out[2]-out2[2]);
                                //a1=nroot(pow(out3[0],2)+pow(out3[1],2),2)*cos(atan2(out3[0],out3[1])-atan2(out2[0]-out[0],out2[1]-out[1]))+nroot(pow(out4[0],2)+pow(out4[1],2),2)*cos(atan2(out4[0],out4[1])-atan2(out[0]-out2[0],out[1]-out2[1]));
                                a1=dot_product(out3,subtract_vectors(out2,out))+dot_product(out4,subtract_vectors(out,out2));
                                if (nroot(pow(out[0]-out2[0],2)+pow(out[1]-out2[1],2)+pow(out[2]-out2[2],2),2.)-2*ball_radius<epsilon && a1>0)
                                {
                                    if(t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
//            if (i==0 && balls[i]._inflight)
//            {
//                std::cout << "Checked balls!" << std::endl;
//            }

            //check if ball hits the baize.
            if (balls[i]._inflight)
            {
                quadratic=qsolve_quadratic(balls[i]._az[2],balls[i]._az[1],balls[i]._az[0]-ball_radius);
                for (int j=0;j<2;j++)
                {
                    if (quadratic[j]==quadratic[j] && quadratic[j]>=epsilon && quadratic[j]<t)
                    {
                        t0=quadratic[j];
                        out=balls[i].get_position(t0);
                        out2=balls[i].get_velocity(t0);
                        if (out[2]<ball_radius && out2[2]<0)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                if (out[2]<ball_radius && out2[2]<0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                if (out[2]<ball_radius && out2[2]<0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
//            if (i==0 && balls[i]._inflight)
//            {
//                std::cout << "Checked bed!" << std::endl;
//            }

            //broad check for cushions - eliminate doomed cases.
            if (xmin>xlim[0] && xmax<xlim[1] && ymin>ylim[0] && ymax<ylim[1])
            {
                continue;
            }

            //check for collision with cushions.

            //straight x.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ax[2],balls[i]._ax[1],balls[i]._ax[0]-xlim[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=epsilon && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out=balls[i].get_position(t0);
                        out2=balls[i].get_velocity(t0);
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //straight y.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ay[2],balls[i]._ay[1],balls[i]._ay[0]-ylim[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=epsilon && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out=balls[i].get_position(t0);
                        out2=balls[i].get_velocity(t0);
                        if (j==0) {a2=-out2[1];}
                        else {a2=out2[1];}
                        if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
//            if (i==0 && balls[i]._inflight)
//            {
//                std::cout << "Checked straight cushions!" << std::endl;
//            }

//            //stricter check!
//            if (!((ymin<=ylim[0] || ymax>=ylim[1]) && ((xmin<=xlim[0] || xmax>=xlim[1]) || (xmin<=blue_x+4.344 && xmax>=blue_x-4.344))))
//            {
//                continue;
//            }

            //check if hits rim.
            if (balls[i]._inflight)
            {
                for (int j=0;j<2;j++)
                {
                    //solve the octic(!)
                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
                    x2[3]=2*balls[i]._ax[2]*balls[i]._ax[1];
                    x2[2]=2*balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
                    x2[1]=2*balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0]);
                    x2[0]=(balls[i]._ax[0]-mpockets[j][0])*(balls[i]._ax[0]-mpockets[j][0]);
                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
                    y2[3]=2*balls[i]._ay[2]*balls[i]._ay[1];
                    y2[2]=2*balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
                    y2[1]=2*balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]);
                    y2[0]=(balls[i]._ay[0]-mpockets[j][1])*(balls[i]._ay[0]-mpockets[j][1]);
                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
                    z2[3]=2*balls[i]._az[2]*balls[i]._az[1];
                    z2[2]=2*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
                    z2[1]=2*balls[i]._az[1]*balls[i]._az[0];
                    z2[0]=balls[i]._az[0]*balls[i]._az[0];

                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
                    toctic[7]=(2*x2[4]*x2[3])+(2*y2[4]*y2[3])+(2*z2[4]*z2[3]);
                    toctic[6]=(2*x2[4]*x2[2]+x2[3]*x2[3])+(2*y2[4]*y2[2]+y2[3]*y2[3])+(2*z2[4]*z2[2]+z2[3]*z2[3]);
                    toctic[5]=(2*x2[4]*x2[1]+2*x2[3]*x2[2])+(2*y2[4]*y2[1]+2*y2[3]*y2[2])+(2*z2[4]*z2[1]+2*z2[3]*z2[2]);
                    toctic[4]=(2*x2[4]*x2[0]+2*x2[3]*x2[1]+x2[2]*x2[2])+(2*y2[4]*y2[0]+2*y2[3]*y2[1]+y2[2]*y2[2])+(2*z2[4]*z2[0]+2*z2[3]*z2[1]+z2[2]*z2[2]);
                    toctic[3]=(2*x2[3]*x2[0]+2*x2[2]*x2[1])+(2*y2[3]*y2[0]+2*y2[2]*y2[1])+(2*z2[3]*z2[0]+2*z2[2]*z2[1]);
                    toctic[2]=(2*x2[2]*x2[0]+x2[1]*x2[1])+(2*y2[2]*y2[0]+y2[1]*y2[1])+(2*z2[2]*z2[0]+z2[1]*z2[1]);
                    toctic[1]=(2*x2[1]*x2[0])+(2*y2[1]*y2[0])+(2*z2[1]*z2[0]);
                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);

                    toctic[8]+=2*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
                    toctic[7]+=2*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
                    toctic[6]+=2*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
                    toctic[5]+=2*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
                    toctic[4]+=2*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
                    toctic[3]+=2*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
                    toctic[2]+=2*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
                    toctic[1]+=2*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
                    toctic[0]+=2*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));

                    toctic[4]+=-2*(pow(mpocket_r,2)+pow(ball_radius,2))*(x2[4]+y2[4])+2*(pow(mpocket_r,2)-pow(ball_radius,2))*z2[4];
                    toctic[3]+=-2*(pow(mpocket_r,2)+pow(ball_radius,2))*(x2[3]+y2[3])+2*(pow(mpocket_r,2)-pow(ball_radius,2))*z2[3];
                    toctic[2]+=-2*(pow(mpocket_r,2)+pow(ball_radius,2))*(x2[2]+y2[2])+2*(pow(mpocket_r,2)-pow(ball_radius,2))*z2[2];
                    toctic[1]+=-2*(pow(mpocket_r,2)+pow(ball_radius,2))*(x2[1]+y2[1])+2*(pow(mpocket_r,2)-pow(ball_radius,2))*z2[1];
                    toctic[0]+=-2*(pow(mpocket_r,2)+pow(ball_radius,2))*(x2[0]+y2[0])+2*(pow(mpocket_r,2)-pow(ball_radius,2))*z2[0];

                    toctic[0]+=pow(pow(mpocket_r,2)-pow(ball_radius,2),2);

                    for (int k=0;k<8;k++)
                    {
                        monoctic[k]=toctic[k]/toctic[8];
                    }
                    octic=qsolve_octic(monoctic);
                    //now to check the roots.
                    for (int k=0;k<8;k++)
                    {
                        if (octic[k]==octic[k] && octic[k]>=epsilon && octic[k]<t)
                        {
                            t0=octic[k];
                            out=balls[i].get_position(t0);
                            out2=balls[i].get_velocity(t0);
                            a1=nroot(pow(out[2],2)+pow(mpocket_r-nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2),2),2);
                            phi=acos(out[2]/a1);
                            nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                            parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                            vtheta=parspeed/nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                            vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                            if (a1<ball_radius && normal>0)
                            {
                                t=t0;
                                collision=true;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                a1=nroot(pow(out[2],2)+pow(mpocket_r-nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2),2),2);
                                phi=acos(out[2]/a1);
                                nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                vtheta=parspeed/nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                                vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                                if (a1<ball_radius && normal>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                a1=nroot(pow(out[2],2)+pow(mpocket_r-nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2),2),2);
                                phi=acos(out[2]/a1);
                                nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                vtheta=parspeed/nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                                vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                                if (a1<ball_radius && normal>0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
                for (int j=0;j<4;j++)
                {
                    //solve the octic(!)
                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
                    x2[3]=2*balls[i]._ax[2]*balls[i]._ax[1];
                    x2[2]=2*balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
                    x2[1]=2*balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0]);
                    x2[0]=(balls[i]._ax[0]-cpockets[j][0])*(balls[i]._ax[0]-cpockets[j][0]);
                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
                    y2[3]=2*balls[i]._ay[2]*balls[i]._ay[1];
                    y2[2]=2*balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
                    y2[1]=2*balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]);
                    y2[0]=(balls[i]._ay[0]-cpockets[j][1])*(balls[i]._ay[0]-cpockets[j][1]);
                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
                    z2[3]=2*balls[i]._az[2]*balls[i]._az[1];
                    z2[2]=2*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
                    z2[1]=2*balls[i]._az[1]*balls[i]._az[0];
                    z2[0]=balls[i]._az[0]*balls[i]._az[0];

                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
                    toctic[7]=(2*x2[4]*x2[3])+(2*y2[4]*y2[3])+(2*z2[4]*z2[3]);
                    toctic[6]=(2*x2[4]*x2[2]+x2[3]*x2[3])+(2*y2[4]*y2[2]+y2[3]*y2[3])+(2*z2[4]*z2[2]+z2[3]*z2[3]);
                    toctic[5]=(2*x2[4]*x2[1]+2*x2[3]*x2[2])+(2*y2[4]*y2[1]+2*y2[3]*y2[2])+(2*z2[4]*z2[1]+2*z2[3]*z2[2]);
                    toctic[4]=(2*x2[4]*x2[0]+2*x2[3]*x2[1]+x2[2]*x2[2])+(2*y2[4]*y2[0]+2*y2[3]*y2[1]+y2[2]*y2[2])+(2*z2[4]*z2[0]+2*z2[3]*z2[1]+z2[2]*z2[2]);
                    toctic[3]=(2*x2[3]*x2[0]+2*x2[2]*x2[1])+(2*y2[3]*y2[0]+2*y2[2]*y2[1])+(2*z2[3]*z2[0]+2*z2[2]*z2[1]);
                    toctic[2]=(2*x2[2]*x2[0]+x2[1]*x2[1])+(2*y2[2]*y2[0]+y2[1]*y2[1])+(2*z2[2]*z2[0]+z2[1]*z2[1]);
                    toctic[1]=(2*x2[1]*x2[0])+(2*y2[1]*y2[0])+(2*z2[1]*z2[0]);
                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);

                    toctic[8]+=2*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
                    toctic[7]+=2*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
                    toctic[6]+=2*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
                    toctic[5]+=2*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
                    toctic[4]+=2*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
                    toctic[3]+=2*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
                    toctic[2]+=2*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
                    toctic[1]+=2*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
                    toctic[0]+=2*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));

                    toctic[4]+=-2*(pow(cpocket_r,2)+pow(ball_radius,2))*(x2[4]+y2[4])+2*(pow(cpocket_r,2)-pow(ball_radius,2))*z2[4];
                    toctic[3]+=-2*(pow(cpocket_r,2)+pow(ball_radius,2))*(x2[3]+y2[3])+2*(pow(cpocket_r,2)-pow(ball_radius,2))*z2[3];
                    toctic[2]+=-2*(pow(cpocket_r,2)+pow(ball_radius,2))*(x2[2]+y2[2])+2*(pow(cpocket_r,2)-pow(ball_radius,2))*z2[2];
                    toctic[1]+=-2*(pow(cpocket_r,2)+pow(ball_radius,2))*(x2[1]+y2[1])+2*(pow(cpocket_r,2)-pow(ball_radius,2))*z2[1];
                    toctic[0]+=-2*(pow(cpocket_r,2)+pow(ball_radius,2))*(x2[0]+y2[0])+2*(pow(cpocket_r,2)-pow(ball_radius,2))*z2[0];

                    toctic[0]+=pow(pow(cpocket_r,2)-pow(ball_radius,2),2);

                    for (int k=0;k<8;k++)
                    {
                        monoctic[k]=toctic[k]/toctic[8];
                    }
                    octic=qsolve_octic(monoctic);
                    //now to check the roots.
                    for (int k=0;k<8;k++)
                    {
                        if (octic[k]==octic[k] && octic[k]>=epsilon && octic[k]<t)
                        {
                            t0=octic[k];
                            out=balls[i].get_position(t0);
                            out2=balls[i].get_velocity(t0);
                            a1=nroot(pow(out[2],2)+pow(cpocket_r-nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2),2),2);
                            phi=acos(out[2]/a1);
                            nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                            parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                            vtheta=parspeed/nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                            vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                            if (a1<ball_radius && normal>0)
                            {
                                t=t0;
                                collision=true;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                a1=nroot(pow(out[2],2)+pow(cpocket_r-nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2),2),2);
                                phi=acos(out[2]/a1);
                                nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                vtheta=parspeed/nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                                vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                                if (a1<ball_radius && normal>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                a1=nroot(pow(out[2],2)+pow(cpocket_r-nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2),2),2);
                                phi=acos(out[2]/a1);
                                nspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                parspeed=nroot(pow(out2[0],2)+pow(out2[1],2),2)*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
                                vtheta=parspeed/nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                                vphi=nroot(pow(nspeed,2)+pow(out2[2],2),2)/a1;
                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2);
                                if (a1<ball_radius && normal>0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                for (int j=0;j<2;j++)
                {
                    //middle pockets.
                    t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                    t3=2*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2*(balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2*(balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-mpockets[j][0],2)+pow(balls[i]._ay[0]-mpockets[j][1],2)-pow(mpocket_r,2);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=epsilon && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out=balls[i].get_position(t0);
                            a1=nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                            a2=nroot(pow(mpocket_r-a1,2)+pow(ball_radius,2),2);
                            if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                t=t0;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                a1=nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                                a2=nroot(pow(mpocket_r-a1,2)+pow(ball_radius,2),2);
                                if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                a1=nroot(pow(out[0]-mpockets[j][0],2)+pow(out[1]-mpockets[j][1],2),2);
                                a2=nroot(pow(mpocket_r-a1,2)+pow(ball_radius,2),2);
                                if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
                for (int j=0;j<4;j++)
                {
                    //corner pockets.
                    t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                    t3=2*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2*(balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2*(balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-cpockets[j][0],2)+pow(balls[i]._ay[0]-cpockets[j][1],2)-pow(cpocket_r,2);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=epsilon && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out=balls[i].get_position(t0);
                            a1=nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                            a2=nroot(pow(cpocket_r-a1,2)+pow(ball_radius,2),2);
                            if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                t=t0;
                            }
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                a1=nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                                a2=nroot(pow(cpocket_r-a1,2)+pow(ball_radius,2),2);
                                if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                a1=nroot(pow(out[0]-cpockets[j][0],2)+pow(out[1]-cpockets[j][1],2),2);
                                a2=nroot(pow(cpocket_r-a1,2)+pow(ball_radius,2),2);
                                if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

//            if (i==0 && balls[i]._inflight)
//            {
//                std::cout << "Checked rim!" << std::endl;
//            }

            //check curved pocket cushions.
            for (int j=0;j<6;j++)
            {
                if (cush[j]._p1==0)
                {
                    ctemp[0]=cp[0];
                    ctemp[1]=cp[1];
                    ctemp[2]=cp[2];
                }
                else
                {
                    ctemp[0]=mp[0];
                    ctemp[1]=mp[1];
                    ctemp[2]=mp[2];
                }
                if (cush[j]._p2==0)
                {
                    ctemp[3]=cp[0];
                    ctemp[4]=cp[1];
                    ctemp[5]=cp[2];
                }
                else
                {
                    ctemp[3]=mp[0];
                    ctemp[4]=mp[1];
                    ctemp[5]=mp[2];
                }
                for (int k=3;k<6;k++)
                {
                    ctemp[k][0]=cush[j]._length-ctemp[k][0];
                    t0=ctemp[k][3];
                    ctemp[k][3]=-ctemp[k][4];
                    ctemp[k][4]=-t0;
                }
                for (int k=0;k<6;k++)
                {
                    x=cush[j]._x+(ctemp[k][0]*cos(cush[j]._angle)+ctemp[k][1]*sin(cush[j]._angle));
                    y=cush[j]._y+(-ctemp[k][0]*sin(cush[j]._angle)+ctemp[k][1]*cos(cush[j]._angle));
                    t4=pow(balls[i]._ax[2],2)+pow(balls[i]._ay[2],2);
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-x)+balls[i]._ay[2]*(balls[i]._ay[0]-y))+pow(balls[i]._ax[1],2)+pow(balls[i]._ay[1],2);
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-x)+balls[i]._ay[1]*(balls[i]._ay[0]-y));
                    t0=pow(balls[i]._ax[0]-x,2)+pow(balls[i]._ay[0]-y,2)-pow(ctemp[k][2],2);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    //check the validity of quartic solutions.
                    for (int l=0;l<4;l++)
                    {
                        if (quartic[l]==quartic[l] && quartic[l]>=epsilon && quartic[l]<t)
                        {
                            //verify the time of collision.
                            t0=quartic[l];
                            out=balls[i].get_position(t0);
                            out2=balls[i].get_velocity(t0);
                            //check the angle with the curve.
                            a1=atan2(out[0]-x,out[1]-y);
                            std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                            a2=-nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(out2[0],out2[1])-b2);
                            a3=ctemp[k][3]+cush[j]._angle;
                            a3=a3-2*pi*floor(a3/(2*pi));
                            b3=ctemp[k][4]+cush[j]._angle;
                            b3=b3-2*pi*floor(b3/(2*pi));
                            if (a3>pi) {a3=a3-2*pi;}
                            if (a3<-pi) {a3=a3+2*pi;}
                            if (b3>pi) {b3=b3-2*pi;}
                            if (b3<-pi) {b3=b3+2*pi;}
                            if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && fabs(b1-ball_radius)<epsilon && a2>0)
                            {
                                t=t0;
                                collision=true;
                            }
                            else
                            {
                                c=0;
                                while (true)
                                {
                                    c+=1;
                                    out=balls[i].get_position(t0-c*epsilon);
                                    out2=balls[i].get_velocity(t0-c*epsilon);
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && fabs(b1-ball_radius)<epsilon && a2>0)
                                    {
                                        if (t0-c*epsilon>=epsilon)
                                        {
                                            t=t0-c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    out=balls[i].get_position(t0+c*epsilon);
                                    out2=balls[i].get_velocity(t0+c*epsilon);
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && fabs(b1-ball_radius)<epsilon && a2>0)
                                    {
                                        if (t0+c*epsilon<t)
                                        {
                                            t=t0+c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    if (c>1000)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

//            if (i==0 && balls[i]._inflight)
//            {
//                std::cout << "Checked curved cushions!" << std::endl;
//            }

            //check middle straight bits.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ax[2],balls[i]._ax[1],balls[i]._ax[0]-mpline[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=epsilon && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out=balls[i].get_position(t0);
                        out2=balls[i].get_velocity(t0);
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //check corner straight bits.
            for (int j=0;j<8;j++)
            {
                a1=0.5*balls[i]._ax[2]-0.5*cpline[j][0]*balls[i]._ay[2];
                a2=0.5*balls[i]._ax[1]-0.5*cpline[j][0]*balls[i]._ay[1];
                a3=0.5*balls[i]._ax[0]-0.5*cpline[j][0]*balls[i]._ay[0]+0.5*cpline[j][0]*cpline[j][1];
                b1=0.5*balls[i]._ay[2]-0.5*cpline[j][0]*balls[i]._ax[2];
                b2=0.5*balls[i]._ay[1]-0.5*cpline[j][0]*balls[i]._ax[1];
                b3=0.5*balls[i]._ay[0]-0.5*cpline[j][0]*balls[i]._ax[0]-0.5*cpline[j][1];

                t4=pow(a1,2)+pow(b1,2);
                t3=2.*(a1*a2+b1*b2);
                t2=2.*(a1*a3+b1*b3)+pow(a2,2)+pow(b2,2);
                t1=2.*(a2*a3+b2*b3);
                t0=pow(a3,2)+pow(b3,2)-pow(ball_radius,2);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=epsilon && quartic[k]<t)
                    {
                        t0=quartic[k];
                        out=balls[i].get_position(t0);
                        out2=balls[i].get_velocity(t0);
                        a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                        b1=cpline[j][0]*a1+cpline[j][1];
                        a2=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                        a3=nroot(pow(out[0]-a1,2)+pow(out[1]-b1,2),2);
                        if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0 && a3<ball_radius)
                        {
                            t=t0;
                            collision=true;
                        }
                        else
                        {
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out=balls[i].get_position(t0-c*epsilon);
                                out2=balls[i].get_velocity(t0-c*epsilon);
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=nroot(pow(out[0]-a1,2)+pow(out[1]-b1,2),2);
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0 && a3<ball_radius)
                                {
                                    if (t0-c*epsilon>=epsilon)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out=balls[i].get_position(t0+c*epsilon);
                                out2=balls[i].get_velocity(t0+c*epsilon);
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=nroot(pow(out2[0],2)+pow(out2[1],2),2)*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=nroot(pow(out[0]-a1,2)+pow(out[1]-b1,2),2);
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0 && a3<ball_radius)
                                {
                                    if (t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                if (c>1000)
                                {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i=0;i<22;i++)
        {
            //deal with balls on the rim.
            if (!balls[i]._onrim || balls[i]._potted)
            {
                continue;
            }

            //check when ball hits cushion rim.
            //work out which pocket ball is rolling in.
            c=0;
            for (int j=0;j<2;j++)
            {
                a1=nroot(pow(balls[i]._x-mpockets[j][0],2)+pow(balls[i]._y-mpockets[j][1],2),2);
                if (a1<mpocket_r+epsilon)
                {
                    //rolling in this pocket!
                    for (int k=0;k<balls[i]._rimpos.size();k++)
                    {
                        if (balls[i]._times[k]>t)
                        {
                            break;
                        }
                        //check for each point here.
                        if (j==1)
                        {
                            a2=nroot(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2),2);
                            b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                            a3=nroot(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2),2);
                            b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                        }
                        else
                        {
                            a2=nroot(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2),2);
                            b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                            a3=nroot(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2),2);
                            b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                        }
                        if (((a2<ball_radius+2.094 && b2>0) || (a3<ball_radius+2.094 && b3>0)) && balls[i]._times[k]<t)
                        {
                            t=balls[i]._times[k];
                            collision=true;
                            break;
                        }
                    }
                    c=j+1;
                    break;
                }
            }
            if (c==0)
            {
                for (int j=0;j<4;j++)
                {
                    a1=nroot(pow(balls[i]._x-cpockets[j][0],2)+pow(balls[i]._y-cpockets[j][1],2),2);
                    if (a1<cpocket_r+epsilon)
                    {
                        //rolling in this pocket!
                        for (int k=0;k<balls[i]._rimpos.size();k++)
                        {
                            if (balls[i]._times[k]>t)
                            {
                                break;
                            }
                            //check for each point here.
                            //check curved bit first.
                            if (j==0)
                            {
                                a2=nroot(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2),2);
                                b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=nroot(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2),2);
                                b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==1)
                            {
                                a2=nroot(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2),2);
                                b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=nroot(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2),2);
                                b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==2)
                            {
                                a2=nroot(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2),2);
                                b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=nroot(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2),2);
                                b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==3)
                            {
                                a2=nroot(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2),2);
                                b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=nroot(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2),2);
                                b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            if (((a2<ball_radius+4.5 && b2>0) || (a3<ball_radius+4.5 && b3>0)) && balls[i]._times[k]<t)
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                            //check straight bits.
                            t0=0.5*cpline[2*j][0]*(balls[i]._rimpos[k][1]+cpline[2*j][0]*balls[i]._rimpos[k][0]-cpline[2*j][1]);
                            t1=cpline[2*j][0]*t0+cpline[2*j][1];
                            a2=nroot(pow(balls[i]._rimpos[k][0]-t0,2)+pow(balls[i]._rimpos[k][1]-t1,2),2);
                            b2=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t0-balls[i]._rimpos[k][0],t1-balls[i]._rimpos[k][1]));

                            t2=0.5*cpline[2*j+1][0]*(balls[i]._rimpos[k][1]+cpline[2*j+1][0]*balls[i]._rimpos[k][0]-cpline[2*j+1][1]);
                            t3=cpline[2*j+1][0]*t2+cpline[2*j+1][1];
                            a3=nroot(pow(balls[i]._rimpos[k][0]-t2,2)+pow(balls[i]._rimpos[k][1]-t3,2),2);
                            b3=nroot(pow(balls[i]._rimvel[k][0],2)+pow(balls[i]._rimvel[k][1],2),2)*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t2-balls[i]._rimpos[k][0],t3-balls[i]._rimpos[k][1]));
                            if (((a2<ball_radius && t0>cpline[2*j][2] && t0<cpline[2*j][3] && b2>0) || (a3<ball_radius && t2>cpline[2*j+1][2] && t2<cpline[2*j+1][3] && b3>0)) && balls[i]._times[k]<t)
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                        }
                        c=j+3;
                        break;
                    }
                }
            }

            //check for collision with other balls.
            for (int j=0;j<22;j++)
            {
                if (i==j || balls[j]._potted)
                {
                    continue;
                }
                if (balls[j]._onrim)
                {
                    //check which pocket.
                    for (int k=0;k<2;k++)
                    {
                        a1=nroot(pow(balls[j]._x-mpockets[k][0],2)+pow(balls[j]._y-mpockets[k][1],2),2);
                        if (a1<mpocket_r)
                        {
                            if (c==k+1)
                            {
                                //same pocket!
                                //check if and when they collide.
                                for (int m=0;m<std::min(balls[i]._times.size(),balls[j]._times.size());m++)
                                {
                                    if (balls[i]._times[m]>t)
                                    {
                                        break;
                                    }
                                    a1=nroot(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2),2);
                                    //a2=nroot(pow(balls[i]._rimvel[m][0],2)+pow(balls[i]._rimvel[m][1],2),2)*cos(atan2(balls[i]._rimvel[m][0],balls[i]._rimvel[m][1])-atan2(balls[j]._rimpos[m][0]-balls[i]._rimpos[m][0],balls[j]._rimpos[m][1]-balls[i]._rimpos[m][1]));
                                    //a2+=nroot(pow(balls[j]._rimvel[m][0],2)+pow(balls[j]._rimvel[m][1],2),2)*cos(atan2(balls[j]._rimvel[m][0],balls[j]._rimvel[m][1])-atan2(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],balls[i]._rimpos[m][0]-balls[j]._rimpos[m][1]));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2*ball_radius && a2>0 && balls[i]._times[m]<t)
                                    {
                                        t=balls[i]._times[m];
                                        collision=true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                    for (int k=0;k<4;k++)
                    {
                        a1=nroot(pow(balls[j]._x-cpockets[k][0],2)+pow(balls[j]._y-cpockets[k][1],2),2);
                        if (a1<cpocket_r)
                        {
                            if (c==k+3)
                            {
                                //same pocket!
                                //check if and when they collide.
                                for (int m=0;m<std::min(balls[i]._times.size(),balls[j]._times.size());m++)
                                {
                                    if (balls[i]._times[m]>t)
                                    {
                                        break;
                                    }
                                    a1=nroot(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2),2);
                                    //a2=nroot(pow(balls[i]._rimvel[m][0],2)+pow(balls[i]._rimvel[m][1],2),2)*cos(atan2(balls[i]._rimvel[m][0],balls[i]._rimvel[m][1])-atan2(balls[j]._rimpos[m][0]-balls[i]._rimpos[m][0],balls[j]._rimpos[m][1]-balls[i]._rimpos[m][1]));
                                    //a2+=nroot(pow(balls[j]._rimvel[m][0],2)+pow(balls[j]._rimvel[m][1],2),2)*cos(atan2(balls[j]._rimvel[m][0],balls[j]._rimvel[m][1])-atan2(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],balls[i]._rimpos[m][0]-balls[j]._rimpos[m][1]));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2*ball_radius && a2>0 && balls[i]._times[m]<t)
                                    {
                                        t=balls[i]._times[m];
                                        collision=true;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                else
                {
                    //check bounds.
                    if (c<3)
                    {
                        //middle pocket.
                    }
                    else
                    {
                        //corner pocket.
                    }
                }
            }
        }

        //add positions to list for each timestep.
        start=ceil(totaltime/timestep);
        for (int i=0;i<int(floor((totaltime+t)/timestep)-start)+1;i++)
        {
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted)
                {
                    temp[3*j]=-100;
                    temp[3*j+1]=-100;
                    temp[3*j+2]=-100;
                    continue;
                }
                if (!balls[j]._onrim)
                {
                    out=balls[j].get_position((start+i)*timestep-totaltime);
                    temp[3*j]=out[0];
                    temp[3*j+1]=out[1];
                    temp[3*j+2]=out[2];
                }
                else
                {
                    t0=(start+i)*timestep-totaltime;
                    if (int(floor(1000*t0))==balls[i]._rimpos.size()-1)
                    {
                        out=balls[i]._rimpos[balls[i]._rimpos.size()-1];
                        temp[3*j]=out[0];
                        temp[3*j+1]=out[1];
                        temp[3*j+2]=out[2];
                    }
                    else
                    {
                        //interpolate.
                        a1=balls[j]._times[int(floor(1000*t0)+1)]-balls[j]._times[int(floor(1000*t0))];
                        out=balls[j]._rimpos[int(floor(1000*t0))];
                        out2=balls[j]._rimpos[int(floor(1000*t0)+1)];
                        temp[3*j]=out[0]+(out2[0]-out[0])*(t0-balls[j]._times[int(floor(1000*t0))])/a1;
                        temp[3*j+1]=out[1]+(out2[1]-out[1])*(t0-balls[j]._times[int(floor(1000*t0))])/a1;
                        temp[3*j+2]=out[2]+(out2[2]-out[2])*(t0-balls[j]._times[int(floor(1000*t0))])/a1;
                    }
                }
            }
            pos.push_back(temp);
        }
        totaltime+=t;

        //fast forward to the next event. new velocity and displacement for the balls.
        //make sure the balls are positioned to collide properly (rounding error?)
        //update spins here.

        grid={};
        grid_index={};
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted)
            {
                continue;
            }
            if (!balls[i]._onrim)
            {
                //position.
                out=balls[i].get_position(t);
                balls[i]._x=out[0];
                balls[i]._y=out[1];
                balls[i]._z=out[2];
                if (out[2]<-ball_radius)
                {
                    balls[i]._potted=true;
                    pos[pos.size()-1][3*i]=-100.;
                    pos[pos.size()-1][3*i+1]=-100.;
                    pos[pos.size()-1][3*i+2]=-100.;
                }
                //velocity.
                if (!balls[i]._potted)
                {
                    out=balls[i].get_velocity(t);
                    if (fabs(out[0])>epsilon) {balls[i]._vx=out[0];}
                    else {balls[i]._vx=0.;}
                    if (fabs(out[1])>epsilon) {balls[i]._vy=out[1];}
                    else {balls[i]._vy=0.;}
                    if (fabs(out[2])>epsilon) {balls[i]._vz=out[2];}
                    else {balls[i]._vz=0.;}
                    //spin.
                    out=balls[i].get_spin(t);
                    if (fabs(out[0])>epsilon) {balls[i]._xspin=out[0];}
                    else {balls[i]._xspin=0.;}
                    if (fabs(out[1])>epsilon) {balls[i]._yspin=out[1];}
                    else {balls[i]._yspin=0.;}
                    if (fabs(out[2])>epsilon) {balls[i]._rspin=out[2];}
                    else {balls[i]._rspin=0.;}
                }
            }
            else
            {
                //on the rim.
                if (int(floor(1000*t))==balls[i]._rimpos.size()-1)
                {
                    out=balls[i]._rimpos[balls[i]._rimpos.size()-1];
                    balls[i]._x=out[0];
                    balls[i]._y=out[1];
                    balls[i]._z=out[2];
                    out=balls[i]._rimvel[balls[i]._rimvel.size()-1];
                    balls[i]._vx=out[0];
                    balls[i]._vy=out[1];
                    balls[i]._vz=out[2];
                }
                else
                {
                    //interpolate.
                    a1=balls[i]._times[int(floor(1000*t)+1)]-balls[i]._times[int(floor(1000*t))];
                    out=balls[i]._rimpos[int(floor(1000*t))];
                    out2=balls[i]._rimpos[int(floor(1000*t)+1)];
                    balls[i]._x=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                    balls[i]._y=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                    balls[i]._z=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                    out=balls[i]._rimvel[int(floor(1000*t))];
                    out2=balls[i]._rimvel[int(floor(1000*t)+1)];
                    balls[i]._vx=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                    balls[i]._vy=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                    balls[i]._vz=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000*t))])/a1;
                }
                balls[i]._xspin=balls[i]._vx/ball_radius;
                balls[i]._yspin=balls[i]._vy/ball_radius;
                //work out how rspin changes.
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0;}
                if (fabs(balls[i]._vz)<epsilon) {balls[i]._vz=0;}
                if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0;}
                if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0;}
            }
            //update the grid positions in case there is a collision.
            gxmin=int(floor((balls[i]._x-balls[i]._r)/(2.*balls[i]._r)));
            gxmax=int(floor((balls[i]._x+balls[i]._r)/(2.*balls[i]._r)));
            gymin=int(floor((balls[i]._y-balls[i]._r)/(2.*balls[i]._r)));
            gymax=int(floor((balls[i]._y+balls[i]._r)/(2.*balls[i]._r)));

            balls[i]._gpos[0][0]=gxmin;
            balls[i]._gpos[0][1]=gymin;
            balls[i]._gpos[1][0]=gxmax;
            balls[i]._gpos[1][1]=gymin;
            balls[i]._gpos[2][0]=gxmax;
            balls[i]._gpos[2][1]=gymax;
            balls[i]._gpos[3][0]=gxmin;
            balls[i]._gpos[3][1]=gymax;
            //update grid position.

            for (int j=0;j<4;j++)
            {
                xpos=balls[i]._gpos[j][0];
                ypos=balls[i]._gpos[j][1];
                grid[xpos][ypos][grid_index[xpos][ypos]]=balls[i]._order;
                grid_index[xpos][ypos]+=1;
            }

            if (i==0)
            {
//                std::cout << "x:" << balls[i]._x << std::endl;
//                std::cout << "y:" << balls[i]._y << std::endl;
//                std::cout << "z:" << balls[i]._z << std::endl;
//                std::cout << "vx:" << balls[i]._vx << std::endl;
//                std::cout << "vy:" << balls[i]._vy << std::endl;
//                std::cout << "vz:" << balls[i]._vz << std::endl;
//                std::cout << "xspin:" << double(balls[i]._xspin) << std::endl;
//                std::cout << "yspin:" << double(balls[i]._yspin) << std::endl;
//                std::cout << "rspin:" << double(balls[i]._rspin) << std::endl;
            }
        }

        //perform collision check.
        if (collision)
        {
            do
            {
                dv=collisions(balls,cush);
                std::cout << "Collision" << std::endl;
                for (int i=0;i<22;i++)
                {
                    if (balls[i]._potted)
                    {
                        continue;
                    }
                    balls[i]._vx+=dv(2*i);
                    balls[i]._vy+=dv(2*i+1);
                }
                //std::cout << "dv: " << dv << std::endl;
            } while (dv!=zero);
            for (int i=0;i<22;i++)
            {
                if (balls[i]._potted)
                {
                    continue;
                }
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0;}
            }
        }
    }
    return pos;
}

#endif // OBJECTS_H_INCLUDED
