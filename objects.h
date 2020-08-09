#ifndef OBJECTS_H_INCLUDED
#define OBJECTS_H_INCLUDED

#include <math.h>
#include <Eigen/Dense>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <set>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random.hpp>
#include "polystuff.h"

using namespace boost::numeric::odeint;

//lengths are in inches.

const int framerate=100;
const double timestep=1./framerate;

const double epsilon=pow(10.,-12.);

const int updates=floor(1./(timestep*framerate));

const double in2m=1./39.37;
const double pi=3.1415926535897932384626433832795028841971;
const double gravity=9.81/in2m;

const double table_length=140.25;
const double table_width=70.125;
const double cush_thickness=1.875;
const double rail_thickness=2.445;
const double bounciness=0.7;
const double mus=0.0075; //0.005-0.015
const double muk=0.15; //0.15-0.4

const double panel_ratio=0.2;
const double raw_width=(2.*rail_thickness+2.*cush_thickness+table_length);
const double raw_height=(2.*rail_thickness+2.*cush_thickness+table_width)*(1.+panel_ratio);
const double dfactor=sf::VideoMode::getDesktopMode().height/raw_height;
const double window_width=(2.*rail_thickness+2.*cush_thickness+table_length)*dfactor;
const double window_height=(2.*rail_thickness+2.*cush_thickness+table_width)*dfactor;
const double panel_height=window_height*panel_ratio;

const double ball_radius=0.02625/in2m;
const double ball_mass=0.141;
const double eball=1.; //coefficient of restitution.

const double mpockets[2][2]={{rail_thickness+cush_thickness+table_width,rail_thickness+2.*cush_thickness+table_width+0.156},
                             {rail_thickness+cush_thickness+table_width,rail_thickness-0.156}};
const double cpockets[4][2]={{rail_thickness,rail_thickness},
                             {rail_thickness,rail_thickness+2.*cush_thickness+table_width},
                             {rail_thickness+2.*cush_thickness+table_length,rail_thickness},
                             {rail_thickness+2.*cush_thickness+table_length,rail_thickness+2.*cush_thickness+table_width}};

const double spot_r=0.25;

const double mpocket_r=2.094;
const double cpocket_r=3.5;

const double pround=3.1;
const double k1=8.595-9.*sin(pi/4.);
const double k2=k1+4.5*sin(pi/4.)-2.445;

const double cueball_break_x=rail_thickness+cush_thickness+table_length-29.;
const double cueball_break_y=rail_thickness+cush_thickness+table_width/2.+5.3;

const double yellow_x=rail_thickness+cush_thickness+table_length-29.;
const double yellow_y=rail_thickness+cush_thickness+table_width/2.+11.687;
const double green_x=yellow_x;
const double green_y=rail_thickness+cush_thickness+table_width/2.-11.687;
const double brown_x=yellow_x;
const double brown_y=rail_thickness+cush_thickness+table_width/2.;
const double blue_x=rail_thickness+cush_thickness+table_width;
const double blue_y=brown_y;
const double pink_x=rail_thickness+cush_thickness+table_width/2.;
const double pink_y=brown_y;
const double black_x=rail_thickness+cush_thickness+12.75;
const double black_y=brown_y;

const double colourpos[6][2]={{yellow_x,yellow_y},{green_x,green_y},{brown_x,brown_y},{blue_x,blue_y},{pink_x,pink_y},{black_x,black_y}};

const sf::Color yellow_col=sf::Color(255,255,0);
const sf::Color green_col=sf::Color(0,150,0);
const sf::Color brown_col=sf::Color(131,87,43);
const sf::Color blue_col=sf::Color(0,0,255);
const sf::Color pink_col=sf::Color(255,105,180);
const sf::Color black_col=sf::Color(0,0,0);

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
const double cush_z1=0.03/in2m;
const double cush_z2=0.039/in2m;

//pocket things.
const double cpocket_angle=atan2(3.33402,1.06503);
const double mpocket_angle=atan2(1.93152,0.652752+0.156);

std::array<std::array<double,5>,3> cp={{{0.5*k1,0.5*k1,sqrt(2.)*0.5*k1-ball_radius,0.75*pi,-0.75*pi},
                                        {6.15,-2.445,4.5+ball_radius,-0.25*pi,-atan(6./11.)},
                                        {5.43,-1.125,3.+ball_radius,-atan(6./11.),0.}}};
std::array<std::array<double,5>,3> mp={{{0.,1.875,pround-ball_radius,pi-asin(1.625/pround),pi},
                                        {3.719,-0.438,2.094+ball_radius,-0.5*pi,-atan(0.625/0.812)},
                                        {4.344,-1.25,3.125+ball_radius,-atan(0.625/0.812),0.}}};

std::array<double,2> xlim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_length-ball_radius};
std::array<double,2> ylim={rail_thickness+cush_thickness+ball_radius,rail_thickness+cush_thickness+table_width-ball_radius};

std::array<double,2> mpline={blue_x-1.625,blue_x+1.625};
std::array<std::array<double,4>,8> cpline={{{1.,-k1,rail_thickness+k1,rail_thickness+k2},
                                            {1.,k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+table_width-k1,rail_thickness,rail_thickness+k2-k1},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+table_width+k1,rail_thickness+k1,rail_thickness+k2},
                                            {-1.,2.*rail_thickness+2.*cush_thickness+k1+table_length,rail_thickness+2.*cush_thickness+table_length-k2+k1,rail_thickness+2.*cush_thickness+table_length},
                                            {-1.,2.*rail_thickness+2.*cush_thickness-k1+table_length,rail_thickness+2.*cush_thickness+table_length-k2,rail_thickness+2.*cush_thickness+table_length-k1},
                                            {1.,-k1-table_width,rail_thickness+2.*cush_thickness+table_length-k2+k1,rail_thickness+2.*cush_thickness+table_length},
                                            {1.,k1-table_width,rail_thickness+2.*cush_thickness+table_length-k2,rail_thickness+2.*cush_thickness+table_length-k1}}};

//matrix.
Eigen::Matrix<double,46,46> M_=Eigen::MatrixXd::Zero(46,46);

std::tuple<double,double> add_vectors(double v1,double a1,double v2,double a2)
{
    double x=v1*sin(a1)+v2*sin(a2);
    double y=v1*cos(a1)+v2*cos(a2);

    return std::make_tuple(sqrt(pow(x,2.)+pow(y,2.)),atan2(x,y));
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

double get_relspeed(std::array<double,3> posa,std::array<double,3> va,std::array<double,3> posb,std::array<double,3> vb)
{
    double mag=sqrt(pow(posb[0]-posa[0],2.)+pow(posb[1]-posa[1],2.)+pow(posb[2]-posa[2],2.));
    std::array<double,3> a;
    a[0]=(posb[0]-posa[0])/mag;
    a[1]=(posb[1]-posa[1])/mag;
    a[2]=(posb[2]-posa[2])/mag;

    double relspeedb=a[0]*(va[0]-vb[0])+a[1]*(va[1]-vb[1])+a[2]*(va[2]-vb[2]);
    return relspeedb;
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
        std::array<double,3> _ax={};
        std::array<double,3> _ay={};
        std::array<double,3> _az={};
        std::array<double,2> _avx={};
        std::array<double,2> _avy={};
        std::array<double,2> _avz={};
        std::array<double,2> _awx={};
        std::array<double,2> _awy={};
        std::array<double,2> _awz={};
        std::vector<std::array<double,3>> _rimpos={};
        std::vector<std::array<double,3>> _rimvel={};
        std::vector<double> _times={};
        double _t=999.;//till next phase change.

        //graphics.
        sf::CircleShape _shape;

        //class functions.
        Ball();
        void update_equation();
};

Ball::Ball()
{
    _shape.setRadius(_r*dfactor);
    _shape.setOrigin(dfactor*_r/2.,dfactor*_r/2.);
    _shape.setPosition(sf::Vector2f(-1000.,-1000.));
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

    if (_x<xlim[0] || _x>xlim[1] || _y<ylim[0] || _y>ylim[1])
    {
        for (int i=0;i<2;i++)
        {
            dist=sqrt(pow(_x-mpockets[i][0],2.)+pow(_y-mpockets[i][1],2.));
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
                dist=sqrt(pow(mpocket_r-dist,2.)+pow(_z,2.));
                if (dist<=ball_radius+epsilon)
                {
                    //on the rim.
                    //check if normal force is zero.
                    phi=acos(_z/dist);
                    nspeed=sqrt(pow(_vx,2.)+pow(_vy,2.))*cos(atan2(mpockets[i][0]-_x,mpockets[i][1]-_y)-atan2(_vx,_vy));
                    parspeed=sqrt(pow(_vx,2.)+pow(_vy,2.))*cos(atan2(mpockets[i][0]-_x,mpockets[i][1]-_y)+0.5*pi-atan2(_vx,_vy));
                    vtheta=parspeed/sqrt(pow(_x-mpockets[i][0],2.)+pow(_y-mpockets[i][1],2.));
                    vphi=sqrt(pow(nspeed,2.)+pow(_vz,2.))/dist;
                    normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
                    //+ve parspeed means going clockwise round pocket.
                    if (normal>0. && vphi>=0.)
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

                        integrate_const(stepper,rimodem,xi,0.,0.3,0.001,push_back_state_and_time(output,_times));
                        size_t ind=0;

                        for (int j=0;j<300;j++)
                        {
                            //get upper limit on the time.
                            ind=j;
                            normal=gravity*cos(output[j][0])-ball_radius*pow(output[j][1],2.)+(mpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][0])*pow(output[j][3],2.);
                            //add on the xyz coords here.
                            relx=(mpocket_r-ball_radius*sin(output[j][0]))*sin(output[j][2]);
                            rely=-(mpocket_r-ball_radius*sin(output[j][0]))*cos(output[j][2]);
                            temp={mpockets[i][0]+relx*cos(_angle)+rely*sin(_angle),mpockets[i][1]-relx*sin(_angle)+rely*cos(_angle),ball_radius*cos(output[j][0])};
                            _rimpos.push_back(temp);
                            std::tie(nspeed,theta)=add_vectors((mpocket_r-ball_radius*sin(output[j][0]))*output[j][3],atan2(mpockets[i][0]-temp[0],mpockets[i][1]-temp[1])+0.5*pi,ball_radius*cos(output[j][0])*output[j][1],atan2(mpockets[i][0]-temp[0],mpockets[i][1]-temp[1]));
                            temp={nspeed*sin(theta),nspeed*cos(theta),-ball_radius*sin(output[j][0])*output[j][1]};
                            _rimvel.push_back(temp);
                            if (normal<0. || (output[j][0]<0. && output[j][1]<0.))
                            {
                                break;
                            }
                        }
                        for (int j=0;j<300-ind;j++)
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
            dist=sqrt(pow(_x-cpockets[i][0],2.)+pow(_y-cpockets[i][1],2.));
            if (_order==1)
            {
                if (dist<2*ball_radius)
                {
                    std::cout << "Distance to rim:" << dist-cpocket_r << std::endl;
                }
            }
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
                dist=sqrt(pow(cpocket_r-dist,2.)+pow(_z,2.));
                if (dist<=ball_radius+epsilon)
                {
                    //on the rim.
                    //check if normal force is zero.
                    phi=acos(_z/dist);
                    nspeed=sqrt(pow(_vx,2.)+pow(_vy,2.))*cos(atan2(cpockets[i][0]-_x,cpockets[i][1]-_y)-atan2(_vx,_vy));
                    parspeed=sqrt(pow(_vx,2.)+pow(_vy,2.))*cos(atan2(cpockets[i][0]-_x,cpockets[i][1]-_y)+0.5*pi-atan2(_vx,_vy));
                    vtheta=parspeed/sqrt(pow(_x-cpockets[i][0],2.)+pow(_y-cpockets[i][1],2.));
                    vphi=sqrt(pow(nspeed,2.)+pow(_vz,2.))/dist;
                    normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
                    //+ve parspeed means going clockwise round pocket.
                    if (normal>0. && vphi>=0.)
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

                        integrate_const(stepper,rimodec,xi,0.,0.3,0.001,push_back_state_and_time(output,_times));
                        size_t ind=0;

                        for (int j=0;j<300;j++)
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
                            if (normal<0. || (output[j][0]<0. && output[j][1]<0.))
                            {
                                break;
                            }
                        }
                        for (int j=0;j<300-ind;j++)
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
    }

    if (!_onrim)
    {
        //check if in ballistic motion.
//        if (_z>ball_radius || _vz>0.)
//        {
//            _inflight=true;
//            _az[2]=-0.5*gravity;
//            _avz[1]=-gravity;
//        }

        //on bed of table and not moving up or down.
        if (!_inflight)
        {
            _vz=0.;
            _avz[0]=0.;
            _z=ball_radius;
            _az[0]=ball_radius;
            _awz[1]=-sgn(_rspin)*mus*gravity;
            double m=sqrt(pow(_vx-_r*_xspin,2.)+pow(_vy-_r*_yspin,2.));
            if ((fabs(_xspin*_r-_vx)<epsilon && fabs(_yspin*_r-_vy)<epsilon))
            {
                //both rolling.
                m=sqrt(pow(_vx,2.)+pow(_vy,2.));
            }

            if (fabs(_vx)<epsilon && fabs(_xspin)<epsilon)
            {
                //stopped in x.
                _ax[2]=0.;
                _avx[1]=0.;
                _awx[1]=0.;
            }
            else if (fabs(_xspin*_r-_vx)<epsilon)
            {
                //rolling in x.
                _vx=_xspin*_r;
                _ax[2]=-0.5*mus*gravity*_vx/m;
                _t=fmin(m/(mus*gravity),_t);
                _avx[1]=-mus*gravity*_vx/m;
                _awx[1]=-mus*gravity*_vx/(_r*m);
            }
            else
            {
                //sliding in x.
                _ax[2]=0.5*muk*gravity*(_xspin*_r-_vx)/m;
                _t=fmin((2./7.)*m/(muk*gravity),_t);
                _avx[1]=(_xspin*_r-_vx)*muk*gravity/m;
                _awx[1]=-(_xspin*_r-_vx)*2.5*muk*gravity/(_r*m);
            }

            if (fabs(_vy)<epsilon && fabs(_yspin)<epsilon)
            {
                //stopped in y.
                _ay[2]=0.;
                _avy[1]=0.;
                _awy[1]=0.;
            }
            else if (fabs(_yspin*_r-_vy)<epsilon)
            {
                //rolling in y.
                _vy=_yspin*_r;
                _ay[2]=-0.5*mus*gravity*_vy/m;
                _t=fmin(m/(mus*gravity),_t);
                _avy[1]=-mus*gravity*_vy/m;
                _awy[1]=-mus*gravity*_vy/(_r*m);
            }
            else
            {
                //sliding in y.
                _ay[2]=0.5*muk*gravity*(_yspin*_r-_vy)/m;
                _t=fmin((2./7.)*m/(muk*gravity),_t);
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
                if (theta<0.)
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
                        if (theta<0. && phi>=epsilon)
                        {
                            _t=phi;
                            break;
                        }
                        phi=roots[i]+c*epsilon;
                        theta=_az[2]*phi*phi+_az[1]*phi+_az[0]+ball_radius;
                        if (theta<0. && phi<_t)
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

//    if (_order==1 && _onrim)
//    {
//        std::cout << "On the rim!" << std::endl;
//        for (int i=0;i<1;i++)
//        {
//            std::cout << _rimpos[i][0] << " : " << _rimpos[i][1] << " : " << _rimpos[i][2] << " : " << _times[i] << std::endl;
//            std::cout << _rimvel[i][0] << " : " << _rimvel[i][1] << " : " << _rimvel[i][2] << std::endl;
//        }
//    }
//    if (_order==1 && _inflight)
//    {
//        std::cout << "In flight!" << std::endl;
//    }
//    if (_order==1 && !_inflight && !_onrim && !_potted)
//    {
//        std::cout << "On bed!" << std::endl;
//    }
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
    double k1=8.595-9*sin(pi/4.);
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
            points[1][1]=0.;
            points[2][0]=k1+delta;
            points[2][1]=0.;
            points[3][0]=_length-1.625-0.438;
            points[3][1]=0.;
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
            points[2][1]=0.;
            points[3][0]=_length-(k1+delta);
            points[3][1]=0.;
            points[4][0]=_length-k1;
            points[4][1]=0.;
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
        points[1][1]=0.;
        points[2][0]=k1+delta;
        points[2][1]=0.;
        points[3][0]=_length-(k1+delta);
        points[3][1]=0.;
        points[4][0]=_length-k1;
        points[4][1]=0.;
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
            angle=0.5*pi*i/50.;
            x=cpocket_r*sin(angle);
            y=cpocket_r*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=0.75*pi+pi*(i-51.)/50.;
            x=0.5*k1+0.5*k1*sqrt(2.)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.)*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
    }
    else
    {
        _pocketshape.setPointCount(102);
        for (int i=0;i<51;i++)
        {
            angle=-0.5*pi+pi*i/50.;
            x=mpocket_r*sin(angle);
            y=-0.156+mpocket_r*cos(angle);
            _pocketshape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi-asin(1.625/pround)+2.*asin(1.625/pround)*(i-51.)/50.;
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
            angle=pi*1.25-0.5*pi*i/50.;
            x=0.5*k1+0.5*k1*sqrt(2.)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.)*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi+0.25*pi*(i-51.)/50.;
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
            angle=pi-asin(1.625/pround)*i/50.;
            x=pround*sin(angle);
            y=1.875+pround*cos(angle);
            _p1shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi-asin(1.625/pround)+asin(1.625/pround)*(i-51.)/50.;
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
            angle=1.25*pi-0.5*pi*i/50.;
            x=_length-0.5*k1+0.5*k1*sqrt(2.)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.)*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=0.75*pi+0.25*pi*(i-51.)/50.;
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
            angle=pi+asin(1.625/pround)-asin(1.625/pround)*i/50.;
            x=_length+pround*sin(angle);
            y=1.875+pround*cos(angle);
            _p2shape.setPoint(i,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
        }
        for (int i=51;i<102;i++)
        {
            angle=pi+asin(1.625/pround)*(i-51.)/50.;
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
            angle=pi*5./4.-i*(pi/2.)/50.;
            x=0.5*k1+0.5*k1*sqrt(2.)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.)*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=pi*7./4.+(pi/4.-atan(6./11.))*i/50.;
            x=6.15+4.5*sin(angle);
            y=-2.445+4.5*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=atan(6./11.)*(i/50.-1.);
            x=5.43+3*sin(angle);
            y=-1.125+3.*cos(angle);
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
            angle=pi*3./2.+atan(0.812/0.625)*i/50.;
            x=3.719+2.094*sin(angle);
            y=-0.438+2.094*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<51;i++)
        {
            angle=pi*3./2.+atan(0.812/0.625)+atan(0.625/0.812)*i/50.;
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
            angle=atan(6./11.)*i/50.;
            x=_length-5.43+3.*sin(angle);
            y=-1.125+3.*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<50;i++)
        {
            angle=atan(6./11.)+(pi/4.-atan(6./11.))*i/50.;
            x=_length-6.15+4.5*sin(angle);
            y=-2.445+4.5*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<1;i++)
        {
            angle=pi*5./4.-i*(pi/2.)/50.;
            x=_length-0.5*k1+0.5*k1*sqrt(2.)*sin(angle);
            y=0.5*k1+0.5*k1*sqrt(2.)*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
    }
    else
    {
        for (int i=0;i<50;i++)
        {
            angle=atan(0.625/0.812)*i/50.;
            x=_length-4.344+3.125*sin(angle);
            y=-1.25+3.125*cos(angle);
            _shape.setPoint(pos,sf::Vector2f((x*cos(_angle)+y*sin(_angle)+_x)*dfactor,window_height-dfactor*(-x*sin(_angle)+y*cos(_angle)+_y)));
            pos+=1;
        }
        for (int i=0;i<51;i++)
        {
            angle=atan(0.625/0.812)+(pi/2.-atan(0.625/0.812))*i/50.;
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
    double relx=(x-_x)*cos(_angle)-(y-_y)*sin(_angle);
    double rely=(x-_x)*sin(_angle)+(y-_y)*cos(_angle);

    double min_dist=999.;
    double normal_angle;

    double a1;
    double a2;

    double dist;
    double angle;

    //check pocket closest to the rel origin.
    //corner pocket.
    if (_p1==0)
    {
        if (height<0.)
        {
            angle=atan2(relx,rely);
            if (angle>=0.25*pi && angle<cpocket_angle)
            {
                dist=sqrt(pow(relx,2.)+pow(rely,2.));
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
            dist=fabs(sqrt(pow(k1*0.5-relx,2.)+pow(k1*0.5-rely,2.))-0.5*k1*sqrt(2.));
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
            dist=sqrt(pow(a1-relx,2.)+pow(a2-rely,2.));
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
            dist=fabs(sqrt(pow(relx-6.15,2.)+pow(rely+2.445,2.))-4.5);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check curvy 2.
        angle=atan2(relx-5.43,rely+1.125);
        if (angle>=-atan(6./11.) && angle<0.)
        {
            dist=fabs(sqrt(pow(relx-5.43,2.)+pow(rely+1.125,2.))-3.);
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
                    normal_angle=0.;
                }
            }
            //check curvy.
            angle=atan2(relx-_length+4.344,rely+1.25);
            if (angle>=0. && angle<atan(0.625/0.812))
            {
                dist=fabs(sqrt(pow(relx-_length+4.344,2.)+pow(rely+1.25,2.))-3.125);
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
                dist=fabs(sqrt(pow(relx-_length+3.719,2.)+pow(rely+0.438,2.))-2.094);
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
            if (angle<=asin(1.625/pround) && angle>0.)
            {
                dist=fabs(sqrt(pow(relx-_length,2.)+pow(-rely+1.875,2.))-pround);
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
                    normal_angle=0.;
                }
            }
            //check curvy.
            angle=atan2(relx-_length+5.43,rely+1.125);
            if (angle>=0. && angle<atan(6./11.))
            {
                dist=fabs(sqrt(pow(relx-_length+5.43,2.)+pow(rely+1.125,2.))-3.);
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
                dist=fabs(sqrt(pow(relx-_length+6.15,2.)+pow(rely+2.445,2.))-4.5);
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
                dist=sqrt(pow(relx-a1,2.)+pow(rely-a2,2.));
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
                dist=fabs(sqrt(pow(k1*0.5-rely,2.)+pow(_length-k1*0.5-relx,2.))-0.5*k1*sqrt(2.));
                if (dist<min_dist)
                {
                    min_dist=dist;
                    normal_angle=angle;
                }
            }
            if (height<0.)
            {
                angle=atan2(relx-_length,rely);
                if (angle>=-cpocket_angle && angle<-0.25*pi)
                {
                    dist=sqrt(pow(relx-_length,2.)+pow(rely,2.));
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
        if (angle<=0. && angle>-asin(1.625/pround))
        {
            dist=fabs(sqrt(pow(-rely+1.875,2.)+pow(relx,2.))-pround);
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
        if (angle>=-pi/2. && angle<-atan(0.625/0.812))
        {
            dist=fabs(sqrt(pow(rely+0.438,2.)+pow(relx-3.719,2.))-2.094);
            if (dist<min_dist)
            {
                min_dist=dist;
                normal_angle=angle;
            }
        }
        //check curvy 2.
        angle=atan2(relx-4.344,rely+1.25);
        if (angle>=-atan(0.625/0.812) && angle<0.)
        {
            dist=fabs(sqrt(pow(relx-4.344,2.)+pow(rely+1.25,2.))-3.125);
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
                normal_angle=0.;
            }
        }
        //check curvy.
        angle=atan2(relx-_length+5.43,rely+1.125);
        if (angle>=0. && angle<atan(6./11.))
        {
            dist=fabs(sqrt(pow(relx-_length+5.43,2.)+pow(rely+1.125,2.))-3.);
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
            dist=fabs(sqrt(pow(relx-_length+6.15,2.)+pow(rely+2.445,2.))-4.5);
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
            dist=sqrt(pow(relx-a1,2.)+pow(rely-a2,2.));
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
            dist=fabs(sqrt(pow(k1*0.5-rely,2.)+pow(_length-k1*0.5-relx,2.))-0.5*k1*sqrt(2.));
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
        double _speed=30.;
        double _mass=0.525;
        double _alpha=0.;
        double _offset=0.;
        double _theta=0.;
        double _eta=0.87;

        double _ballv;
        double _ballvz;
        double _ballparspin;
        double _ballperspin;
        double _ballrspin;

        boost::random_device rd;
        boost::mt19937 _gen{boost::mt19937(rd())};
        boost::normal_distribution<> _ndangle{boost::normal_distribution<>(0.0,pi/3600.)};
        boost::normal_distribution<> _ndalpha{boost::normal_distribution<>(0.0,0.003)};
        boost::normal_distribution<> _ndspeed{boost::normal_distribution<>(0.0,0.0085)};
        boost::normal_distribution<> _ndoffset{boost::normal_distribution<>(0.0,0.04)};

        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varangle{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndangle)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varalpha{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndalpha)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varspeed{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndspeed)};
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > _varoffset{boost::variate_generator<boost::mt19937&,boost::normal_distribution<> >(_gen,_ndoffset)};

        sf::Texture _texture;
        sf::Sprite _sprite;

        Cue();
        void shot();
        void perturb();
};

Cue::Cue()
{
    if (!_texture.loadFromFile("cue3.png"))
    {
        std::cout << "Error loading cue texture!" << std::endl;
    }
    _texture.setSmooth(true);
    _sprite.setTexture(_texture);
//    _sprite.scale(58.*dfactor/1369.,58.*dfactor/1369.);
//    _sprite.setOrigin(0.,26.*58.*dfactor/1369.);
    _sprite.scale(57.*dfactor/5213.,57.*dfactor/5213.);
    sf::FloatRect bounds=_sprite.getLocalBounds();
    _sprite.setOrigin(bounds.left,bounds.top+bounds.height*0.5);
}

void Cue::shot()
{
    double dy=_offset*cos(_theta);
    double dx=_offset*sin(_theta);
    double a=(1.+(ball_mass/_mass)+(2.5/pow(ball_radius,2.))*(pow(_offset,2.)));
    double v=(_speed/a)*(1+sqrt(_eta-(1-_eta)*(_mass/ball_mass)*(1+(2.5/pow(ball_radius,2.))*(pow(_offset,2.)))));
    double w=2.5*v*_offset/pow(ball_radius,2.);

    _ballparspin=w*cos(_theta);
    _ballrspin=-w*sin(_theta)*cos(_alpha);
    _ballperspin=w*sin(_theta)*sin(_alpha);

    _ballv=v*cos(_alpha);
    _ballvz=-v*sin(_alpha);
}

void Cue::perturb()
{
    _angle+=_varangle();
    _alpha+=_varalpha();
    _speed+=_varspeed();

    if (_alpha<0.) {_alpha=0.;}
    if (_alpha>0.5*pi) {_alpha=0.5*pi;}

    if (_speed<0) {_speed=0.;}
    if (_speed>100.5*1.2) {_speed=100.5*1.2;}

    double dx=_offset*sin(_theta);
    double dy=_offset*cos(_theta);

    dx+=_varoffset();
    dy+=_varoffset();

    _offset=sqrt(pow(dx,2.)+pow(dy,2.));
    _theta=atan2(dx,dy);

    if (_offset>ball_radius) {_offset=ball_radius;}
}

class Server
{
    public:
        bool running=true;

        bool player_turn=0;

        std::string names[2]={"PLAYER 1","PLAYER 2"};
        int scores[2]={0,0};
        int frames[2]={0,0};
        int framesbestof=35;
        bool placing_white=true;
        bool touching=false;
        bool isyourturn=true;
        bool isfoul=false;
        bool ismiss=false;
        bool ispush=false;
        bool isredon=true;
        int foulscore=0;
        int nom_colour_order=0; //1 yellow,2 green, etc.
        bool isfreeball=false;

        bool gameover=false;
        int highbreak[2]={0,0};
        int centuries[2]={0,0};
        int current_break=0;

        //info for each shot.
        std::vector<int> ball_hit_order;
        std::vector<int> ball_potted_order;

        sf::IpAddress serverIp;
        unsigned short port;
        sf::TcpListener listener;

        std::array<sf::TcpSocket,2> players;
        std::array<sf::TcpSocket,4> spectators;

        sf::Packet packet;

        //for objects just need balls and cushions.

        //prepare cushions.
        Cushion servercushions[6];

        Ball serverballs[22];

        //grid squares.
        std::array<std::array<std::array<int,7>,39>,73> grid={};
        std::array<std::array<int,39>,73> grid_index={};

        //matrix.
        Eigen::Matrix<double,46,46> sM_=Eigen::MatrixXd::Zero(46,46);

        std::vector<std::array<double,66> > result;

        Server();
        std::vector<std::array<double,66> > simulate(Ball balls[22],Cushion cush[6]);
        Eigen::Matrix<double,46,1> collisions(Ball b[22], Cushion cush[6]);
        void handleIncomingConnections();
        void executionThread();
        void shutdown();
        void respot();
        void rackballs();
        void turnpacket();
        void resetframe();
        bool is_snookered();
};

Server::Server()
{
    for (int i=0;i<2;i++)
    {
        players[i].setBlocking(false);
    }
    for (int i=0;i<4;i++)
    {
        spectators[i].setBlocking(false);
    }

    servercushions[0]=Cushion(mpockets[0][0],mpockets[0][1]-0.156,pi,1,0);
    servercushions[1]=Cushion(mpockets[1][0],mpockets[1][1]+0.156,0.0,1,0);
    servercushions[2]=Cushion(cpockets[0][0],cpockets[0][1],0.0,0,1);
    servercushions[3]=Cushion(cpockets[1][0],cpockets[1][1],pi/2,0,0);
    servercushions[4]=Cushion(cpockets[2][0],cpockets[2][1],pi*3/2,0,0);
    servercushions[5]=Cushion(cpockets[3][0],cpockets[3][1],pi,0,1);

    //prepare matrix.
    for (int i=0;i<44;i++)
    {
        sM_(i,i)=1/ball_mass;
    }
    sM_(44,44)=0;
    sM_(45,45)=0;
}

void Server::resetframe()
{

}

bool Server::is_snookered()
{
    return false;
}

void Server::turnpacket()
{
    //send turn packet.
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            packet << sf::Uint16(2);
            packet << (player_turn==i);
            packet << isfoul << ismiss << placing_white << isredon << isfreeball << gameover;
            packet << sf::Uint32(scores[i]) << sf::Uint32(scores[(i+1)%2]) << sf::Uint32(frames[i]) << sf::Uint32(frames[(i+1)%2]) << sf::Uint32(highbreak[i]) << sf::Uint32(highbreak[(i+1)%2]) << sf::Uint32(centuries[i]) << sf::Uint32(centuries[(i+1)%2]);
            players[i].send(packet);
        }
    }
    for (int i=0;i<4;i++)
    {
        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            packet << sf::Uint16(2);
            packet << false;
            packet << isfoul << ismiss << placing_white << isredon << isfreeball << gameover;
            packet << sf::Uint32(scores[0]) << sf::Uint32(scores[1]) << sf::Uint32(frames[0]) << sf::Uint32(frames[1]) << sf::Uint32(highbreak[0]) << sf::Uint32(highbreak[1]) << sf::Uint32(centuries[0]) << sf::Uint32(centuries[1]);
            spectators[i].send(packet);
        }
    }
}

void Server::rackballs()
{
    //cueball.
    serverballs[0]._x=cueball_break_x;
    serverballs[0]._y=cueball_break_y;
    serverballs[0]._order=1;
    //yellow.
    serverballs[1]._x=yellow_x;
    serverballs[1]._y=yellow_y;
    serverballs[1]._order=2;
    //green.
    serverballs[2]._x=green_x;
    serverballs[2]._y=green_y;
    serverballs[2]._order=3;
    //brown.
    serverballs[3]._x=brown_x;
    serverballs[3]._y=brown_y;
    serverballs[3]._order=4;
    //blue.
    serverballs[4]._x=blue_x;
    serverballs[4]._y=blue_y;
    serverballs[4]._order=5;
    //pink.
    serverballs[5]._x=pink_x;
    serverballs[5]._y=pink_y;
    serverballs[5]._order=6;
    //black.
    serverballs[6]._x=black_x;
    serverballs[6]._y=black_y;
    serverballs[6]._order=7;

    double x;
    double y;
    int i=7;
    int gxmin;
    int gxmax;
    int gymin;
    int gymax;

    for (int row=0;row<5;row++)
    {
        for (int h=-row;h<row+1;h=h+2)
        {
            serverballs[i]._x=pink_x-0.1-2*ball_radius-row*sqrt(3.0)*(ball_radius+DOUBLE_EPSILON);
            serverballs[i]._y=pink_y+h*(ball_radius+DOUBLE_EPSILON);
            serverballs[i]._order=i+1;
            i+=1;
        }
    }

    for (int i=0;i<22;i++)
    {
        serverballs[i]._potted=false;
        serverballs[i]._z=ball_radius;
        serverballs[i]._vx=0.;
        serverballs[i]._vy=0.;
        serverballs[i]._vz=0.;
        serverballs[i]._xspin=0.;
        serverballs[i]._yspin=0.;
        serverballs[i]._rspin=0.;
        gxmin=int(floor((serverballs[i]._x-serverballs[i]._r)/(2.*serverballs[i]._r)));
        gxmax=int(floor((serverballs[i]._x+serverballs[i]._r)/(2.*serverballs[i]._r)));
        gymin=int(floor((serverballs[i]._y-serverballs[i]._r)/(2.*serverballs[i]._r)));
        gymax=int(floor((serverballs[i]._y+serverballs[i]._r)/(2.*serverballs[i]._r)));

        serverballs[i]._gpos[0][0]=gxmin;
        serverballs[i]._gpos[0][1]=gymin;
        serverballs[i]._gpos[1][0]=gxmax;
        serverballs[i]._gpos[1][1]=gymin;
        serverballs[i]._gpos[2][0]=gxmax;
        serverballs[i]._gpos[2][1]=gymax;
        serverballs[i]._gpos[3][0]=gxmin;
        serverballs[i]._gpos[3][1]=gymax;
    }
}

void Server::respot()
{
    double xspot=0.;
    double yspot=0.;
    std::array<bool,6> covered;
    double dx=0.001/in2m;

    for (int i=0;i<6;i++)
    {
        covered[i]=false;
        xspot=colourpos[i][0];
        yspot=colourpos[i][1];
        serverballs[i+1]._z=ball_radius;

        for (int j=0;j<22;j++)
        {
            if (sqrt(pow(xspot-serverballs[j]._x,2.)+pow(yspot-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
            {
                covered[i]=true;
                break;
            }
        }
    }

    for (int i=6;i>0;i--)
    {
        if (!serverballs[i]._potted) {continue;}
        if (!covered[i-1])
        {
            serverballs[i]._x=colourpos[i-1][0];
            serverballs[i]._y=colourpos[i-1][1];
            covered[i-1]=true;
            serverballs[i]._potted=false;
            continue;
        }
        //default spot covered.
        //find highest empty spot.
        for (int j=5;j>-1;j--)
        {
            if (!covered[j])
            {
                serverballs[i]._x=colourpos[j][0];
                serverballs[i]._y=colourpos[j][1];
                covered[j]=true;
                serverballs[i]._potted=false;
                break;
            }
        }
        if (!serverballs[i]._potted) {continue;}

        //move closest to its spot.
        serverballs[i]._x=colourpos[i-1][0];
        serverballs[i]._y=colourpos[i-1][1];
        bool good=true;
        while (true)
        {
            serverballs[i]._x-=dx;
            if (serverballs[i]._x-rail_thickness-cush_thickness<epsilon)
            {
                serverballs[i]._x=-100.;
                break;
            }
            good=true;
            for (int j=0;j<22;j++)
            {
                if (i==j) {continue;}
                if (sqrt(pow(serverballs[i]._x-serverballs[j]._x,2.)+pow(serverballs[i]._y-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
                {
                    good=false;
                    break;
                }
            }
            if (good) {break;}
        }
        if (!serverballs[i]._potted) {continue;}

        serverballs[i]._x=colourpos[i-1][0];
        while (true)
        {
            serverballs[i]._x+=dx;
            good=true;
            for (int j=0;j<22;j++)
            {
                if (i==j) {continue;}
                if (sqrt(pow(serverballs[i]._x-serverballs[j]._x,2.)+pow(serverballs[i]._y-serverballs[j]._y,2.))-2.*ball_radius<epsilon)
                {
                    good=false;
                    break;
                }
            }
            if (good) {break;}
        }
    }
}

Eigen::Matrix<double,46,1> Server::collisions(Ball b[22],Cushion cush[6])
{
    //broad phase collision check.
    Eigen::MatrixXd dv=Eigen::MatrixXd::Zero(46,1);

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
            if (b[i]._vx!=0. || b[i]._vy!=0. || b[i]._vz!=0.)
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
                                    dist=sqrt(pow(dx,2.)+pow(dy,2.));
                                    relspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(dx,dy)-atan2(b[i]._vx,b[i]._vy))-sqrt(pow(b[thing-1]._vx,2.)+pow(b[thing-1]._vy,2.))*cos(atan2(dx,dy)-atan2(b[thing-1]._vx,b[thing-1]._vy));

                                    if (dist-2.*ball_radius<epsilon && relspeed>0.)
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
                                        col+=1;
                                    }
                                }
                            }
                        }
                    }
                }
                //check cushion collisions.
                if (fabs(blue_x-b[i]._x)>table_width-ball_radius-2.*epsilon || fabs(blue_y-b[i]._y)>0.5*table_width-ball_radius-2.*epsilon)
                {
                    double min_dist=999.;
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
                    relspeed=-sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(normal_angle-atan2(b[i]._vx,b[i]._vy));
                    if (min_dist-ball_radius<epsilon && relspeed>0.)
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
                        if (relspin>0.)
                        {
                            //topspin.
                            dw=-5.*relspeed*(muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin<0.)
                        {
                            //backspin.
                            dw=-5.*relspeed*(-muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin==0.)
                        {
                            dw=-5.*relspeed*cushion_diff/pow(ball_radius,2.);
                        }
                        dw=fmax(dw,-relspin-relspeed/ball_radius);
                        b[i]._xspin+=dw*sin(normal_angle+pi);
                        b[i]._yspin+=dw*cos(normal_angle+pi);
                        //adjust sidespin off cushion.
                        double parspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*sin(atan2(b[i]._vx,b[i]._vy)-normal_angle);
                        double dvpar;

                        if ((parspeed>=0. && b[i]._rspin*ball_radius>parspeed) || (parspeed<0. && b[i]._rspin*ball_radius>parspeed))
                        {
                            dw=-5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius<parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        else if ((parspeed>=0. && b[i]._rspin*ball_radius<parspeed) || (parspeed<0. && b[i]._rspin*ball_radius<parspeed))
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
                        std::tie(_speed,_angle)=add_vectors(sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.)),atan2(b[i]._vx,b[i]._vy),dvpar,normal_angle+0.5*pi);
                        b[i]._vx=_speed*sin(_angle);
                        b[i]._vy=_speed*cos(_angle);
                    }
                }
            }
        }
    }

    //simultaneous.
    if (col>0)
    {
        Eigen::MatrixXd v(46,1);
        Eigen::MatrixXd J(col,1);

        for (int i=0;i<22;i++)
        {
            v(2*i,0)=b[i]._vx;
            v(2*i+1,0)=b[i]._vy;
        }
        v(44,0)=0;
        v(45,0)=0;

        //J=(Z.transpose()*M_*Z).colPivHouseholderQr().solve(-(1.+eball)*Z.transpose()*v);
        J=(Z.transpose()*sM_*Z).householderQr().solve(-(1.+eball)*Z.transpose()*v);
        dv=sM_*Z*J;
    }
    return dv;
}

std::vector<std::array<double,66> > Server::simulate(Ball balls[22],Cushion cush[6])
{
    //assumes that the input is not a static scenario (where all balls are still).
    std::vector<std::array<double,66> > pos;
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
        dv(i,0)=0.;
        zero(i,0)=0.;
    }

    //set up the equations for each ball.

    double t; //time until next event.
    double totaltime=0.; //total time elapsed.
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

    bool collision=false;

    ball_hit_order.clear();
    ball_potted_order.clear();
    bool hitdone=false;

    while (true)
    {
//        std::cout << "Time: " << totaltime << std::endl;

//        if (totaltime>0.2)
//        {
//            break;
//        }

        collision=false;
        t=999.;
        for (int i=0;i<22;i++)
        {
            if (balls[i]._potted==false)
            {
                balls[i].update_equation();
                t=fmin(t,balls[i]._t);
            }
        }
        if (t>998.)
        {
            //add final positions.
            for (int i=0;i<22;i++)
            {
                if (!balls[i]._potted)
                {
                    temp[3*i]=balls[i]._x;
                    temp[3*i+1]=balls[i]._y;
                    temp[3*i+2]=balls[i]._z;
                }
                else
                {
                    temp[3*i]=-100.;
                    temp[3*i+1]=-100.;
                    temp[3*i+2]=-100.;
                }
            }
            pos.push_back(temp);
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
            if (balls[i]._vx==0. && balls[i]._vy==0. && balls[i]._vz==0. && balls[i]._xspin==0. && balls[i]._yspin==0.)
            {
                continue;
            }

            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmin-=1.01*ball_radius;
            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmax+=1.01*ball_radius;
            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymin-=1.01*ball_radius;
            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymax+=1.01*ball_radius;

            //check for ball-ball collisions.
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted || balls[j]._onrim || i==j)
                {
                    continue;
                }

                //initial check for i.
                xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmin2-=1.01*ball_radius;
                xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmax2+=1.01*ball_radius;
                ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymin2-=1.01*ball_radius;
                ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymax2+=1.01*ball_radius;

                if (ymin>ymax2 || ymin2>ymax || xmin>xmax2 || xmin2>xmax)
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

                t0=pow(c1,2.)+pow(c2,2.)+pow(c3,2.)-4.*pow(ball_radius,2.);
                t1=2.*(b1*c1+b2*c2+b3*c3);
                t2=2.*(a1*c1+a2*c2+a3*c3)+pow(b1,2.)+pow(b2,2.)+pow(b3,2.);
                t3=2.*(a1*b1+a2*b2+a3*b3);
                t4=pow(a1,2.)+pow(a2,2.)+pow(a3,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        //verify the time of actual collision to be certain.
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                        out2[0]=balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0];
                        out2[1]=balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0];
                        out2[2]=balls[j]._az[2]*t0*t0+balls[j]._az[1]*t0+balls[j]._az[0];
                        out3[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        out3[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        out3[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
                        out4[0]=balls[j]._avx[1]*t0+balls[j]._avx[0];
                        out4[1]=balls[j]._avy[1]*t0+balls[j]._avy[0];
                        out4[2]=balls[j]._avz[1]*t0+balls[j]._avz[0];
                        a1=get_relspeed(out,out3,out2,out4);
                        if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                        {
                            //collision at the specified time.
                            t=t0;
                            collision=true;

                            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmin-=1.01*ball_radius;
                            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmax+=1.01*ball_radius;
                            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymin-=1.01*ball_radius;
                            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymax+=1.01*ball_radius;
                        }
                        else
                        {
                            //adjust the time minutely to ensure a collision.
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ax[1]*(t0-c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ay[1]*(t0-c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._az[1]*(t0-c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0-c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0-c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0-c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ax[1]*(t0+c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ay[1]*(t0+c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._az[1]*(t0+c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0+c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0+c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0+c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if(t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
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

//            //check if ball hits the baize.
//            if (balls[i]._inflight)
//            {
//                quadratic=qsolve_quadratic(balls[i]._az[2],balls[i]._az[1],balls[i]._az[0]-ball_radius);
//                for (int j=0;j<2;j++)
//                {
//                    if (quadratic[j]==quadratic[j] && quadratic[j]>=0. && quadratic[j]<t)
//                    {
//                        t0=quadratic[j];
//                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                        out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                        out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                        //check if not in a pocket.
//                        a1=999.;
//                        a2=999.;
//                        for (int k=0;k<4;k++)
//                        {
//                            //corner check
//                            a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                        }
//                        for (int k=0;k<2;k++)
//                        {
//                            a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                        }
//                        if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                        {
//                            t=t0;
//                            collision=true;
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                //check if not in a pocket.
//                                a1=999.;
//                                a2=999.;
//                                for (int k=0;k<4;k++)
//                                {
//                                    //corner check
//                                    a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                                }
//                                for (int k=0;k<2;k++)
//                                {
//                                    a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                                }
//                                if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                //check if not in a pocket.
//                                a1=999.;
//                                a2=999.;
//                                for (int k=0;k<4;k++)
//                                {
//                                    //corner check
//                                    a1=fmin(a1,sqrt(pow(out[0]-cpockets[k][0],2.)+pow(out[1]-cpockets[k][1],2.)));
//                                }
//                                for (int k=0;k<2;k++)
//                                {
//                                    a2=fmin(a2,sqrt(pow(out[0]-mpockets[k][0],2.)+pow(out[1]-mpockets[k][1],2.)));
//                                }
//                                if (!(a1<cpocket_r+epsilon || a2<mpocket_r+epsilon) && out[2]<ball_radius && out2[2]<0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
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
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
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
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        if (j==0) {a2=-out2[1];}
                        else {a2=out2[1];}
                        if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
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

            //check if hits rim.
            if (balls[i]._inflight)
            {
//                for (int j=0;j<2;j++)
//                {
//                    //solve the octic(!)
//                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
//                    x2[3]=2.*balls[i]._ax[2]*balls[i]._ax[1];
//                    x2[2]=2.*balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
//                    x2[1]=2.*balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0]);
//                    x2[0]=(balls[i]._ax[0]-mpockets[j][0])*(balls[i]._ax[0]-mpockets[j][0]);
//                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
//                    y2[3]=2.*balls[i]._ay[2]*balls[i]._ay[1];
//                    y2[2]=2.*balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
//                    y2[1]=2.*balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]);
//                    y2[0]=(balls[i]._ay[0]-mpockets[j][1])*(balls[i]._ay[0]-mpockets[j][1]);
//                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
//                    z2[3]=2.*balls[i]._az[2]*balls[i]._az[1];
//                    z2[2]=2.*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
//                    z2[1]=2.*balls[i]._az[1]*balls[i]._az[0];
//                    z2[0]=balls[i]._az[0]*balls[i]._az[0];
//
//                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
//                    toctic[7]=(2.*x2[4]*x2[3])+(2.*y2[4]*y2[3])+(2.*z2[4]*z2[3]);
//                    toctic[6]=(2.*x2[4]*x2[2]+x2[3]*x2[3])+(2.*y2[4]*y2[2]+y2[3]*y2[3])+(2.*z2[4]*z2[2]+z2[3]*z2[3]);
//                    toctic[5]=(2.*x2[4]*x2[1]+2.*x2[3]*x2[2])+(2.*y2[4]*y2[1]+2*y2[3]*y2[2])+(2.*z2[4]*z2[1]+2.*z2[3]*z2[2]);
//                    toctic[4]=(2.*x2[4]*x2[0]+2.*x2[3]*x2[1]+x2[2]*x2[2])+(2.*y2[4]*y2[0]+2.*y2[3]*y2[1]+y2[2]*y2[2])+(2.*z2[4]*z2[0]+2.*z2[3]*z2[1]+z2[2]*z2[2]);
//                    toctic[3]=(2.*x2[3]*x2[0]+2.*x2[2]*x2[1])+(2.*y2[3]*y2[0]+2.*y2[2]*y2[1])+(2.*z2[3]*z2[0]+2.*z2[2]*z2[1]);
//                    toctic[2]=(2.*x2[2]*x2[0]+x2[1]*x2[1])+(2.*y2[2]*y2[0]+y2[1]*y2[1])+(2.*z2[2]*z2[0]+z2[1]*z2[1]);
//                    toctic[1]=(2.*x2[1]*x2[0])+(2.*y2[1]*y2[0])+(2.*z2[1]*z2[0]);
//                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);
//
//                    toctic[8]+=2.*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
//                    toctic[7]+=2.*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
//                    toctic[6]+=2.*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
//                    toctic[5]+=2.*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
//                    toctic[4]+=2.*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
//                    toctic[3]+=2.*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
//                    toctic[2]+=2.*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
//                    toctic[1]+=2.*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
//                    toctic[0]+=2.*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));
//
//                    toctic[4]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[4]+y2[4])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[4];
//                    toctic[3]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[3]+y2[3])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[3];
//                    toctic[2]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[2]+y2[2])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[2];
//                    toctic[1]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[1]+y2[1])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[1];
//                    toctic[0]+=-2.*(pow(mpocket_r,2.)+pow(ball_radius,2.))*(x2[0]+y2[0])+2.*(pow(mpocket_r,2.)-pow(ball_radius,2.))*z2[0];
//
//                    toctic[0]+=pow(pow(mpocket_r,2.)-pow(ball_radius,2.),2.);
//
//                    for (int k=0;k<8;k++)
//                    {
//                        monoctic[k]=toctic[k]/toctic[8];
//                    }
//                    octic=qsolve_octic(monoctic);
//                    //now to check the roots.
//                    for (int k=0;k<8;k++)
//                    {
//                        if (octic[k]==octic[k] && octic[k]>=0. && octic[k]<t)
//                        {
//                            t0=octic[k];
//                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
//                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
//                            out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                            a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                            phi=acos(out[2]/a1);
//                            nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                            vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                            if (a1<ball_radius && normal>0.)
//                            {
//                                t=t0;
//                                collision=true;
//                            }
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(mpocket_r-sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(mpockets[j][0]-out[0],mpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(mpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//                for (int j=0;j<4;j++)
//                {
//                    //solve the octic(!)
//                    x2[4]=balls[i]._ax[2]*balls[i]._ax[2];
//                    x2[3]=2.*balls[i]._ax[2]*balls[i]._ax[1];
//                    x2[2]=2.*balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ax[1]*balls[i]._ax[1];
//                    x2[1]=2.*balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0]);
//                    x2[0]=(balls[i]._ax[0]-cpockets[j][0])*(balls[i]._ax[0]-cpockets[j][0]);
//                    y2[4]=balls[i]._ay[2]*balls[i]._ay[2];
//                    y2[3]=2.*balls[i]._ay[2]*balls[i]._ay[1];
//                    y2[2]=2.*balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1])+balls[i]._ay[1]*balls[i]._ay[1];
//                    y2[1]=2.*balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]);
//                    y2[0]=(balls[i]._ay[0]-cpockets[j][1])*(balls[i]._ay[0]-cpockets[j][1]);
//                    z2[4]=balls[i]._az[2]*balls[i]._az[2];
//                    z2[3]=2.*balls[i]._az[2]*balls[i]._az[1];
//                    z2[2]=2.*balls[i]._az[2]*balls[i]._az[0]+balls[i]._az[1]*balls[i]._az[1];
//                    z2[1]=2.*balls[i]._az[1]*balls[i]._az[0];
//                    z2[0]=balls[i]._az[0]*balls[i]._az[0];
//
//                    toctic[8]=(x2[4]*x2[4])+(y2[4]*y2[4])+(z2[4]*z2[4]);
//                    toctic[7]=(2.*x2[4]*x2[3])+(2.*y2[4]*y2[3])+(2.*z2[4]*z2[3]);
//                    toctic[6]=(2.*x2[4]*x2[2]+x2[3]*x2[3])+(2.*y2[4]*y2[2]+y2[3]*y2[3])+(2.*z2[4]*z2[2]+z2[3]*z2[3]);
//                    toctic[5]=(2.*x2[4]*x2[1]+2.*x2[3]*x2[2])+(2.*y2[4]*y2[1]+2.*y2[3]*y2[2])+(2.*z2[4]*z2[1]+2.*z2[3]*z2[2]);
//                    toctic[4]=(2.*x2[4]*x2[0]+2.*x2[3]*x2[1]+x2[2]*x2[2])+(2.*y2[4]*y2[0]+2.*y2[3]*y2[1]+y2[2]*y2[2])+(2.*z2[4]*z2[0]+2.*z2[3]*z2[1]+z2[2]*z2[2]);
//                    toctic[3]=(2.*x2[3]*x2[0]+2.*x2[2]*x2[1])+(2.*y2[3]*y2[0]+2.*y2[2]*y2[1])+(2.*z2[3]*z2[0]+2.*z2[2]*z2[1]);
//                    toctic[2]=(2.*x2[2]*x2[0]+x2[1]*x2[1])+(2.*y2[2]*y2[0]+y2[1]*y2[1])+(2.*z2[2]*z2[0]+z2[1]*z2[1]);
//                    toctic[1]=(2.*x2[1]*x2[0])+(2.*y2[1]*y2[0])+(2.*z2[1]*z2[0]);
//                    toctic[0]=(x2[0]*x2[0])+(y2[0]*y2[0])+(z2[0]*z2[0]);
//
//                    toctic[8]+=2.*((x2[4]*y2[4])+(y2[4]*z2[4])+(z2[4]*x2[4]));
//                    toctic[7]+=2.*((x2[4]*y2[3]+x2[3]*y2[4])+(y2[4]*z2[3]+y2[3]*z2[4])+(z2[4]*x2[3]+z2[3]*x2[4]));
//                    toctic[6]+=2.*((x2[4]*y2[2]+x2[2]*y2[4]+x2[3]*y2[3])+(y2[4]*z2[2]+y2[2]*z2[4]+y2[3]*z2[3])+(z2[4]*x2[2]+z2[2]*x2[4]+z2[3]*x2[3]));
//                    toctic[5]+=2.*((x2[4]*y2[1]+x2[1]*y2[4]+x2[3]*y2[2]+x2[2]*y2[3])+(y2[4]*z2[1]+y2[1]*z2[4]+y2[3]*z2[2]+y2[2]*z2[3])+(z2[4]*x2[1]+z2[1]*x2[4]+z2[3]*x2[2]+z2[2]*x2[3]));
//                    toctic[4]+=2.*((x2[4]*y2[0]+x2[0]*y2[4]+x2[3]*y2[1]+x2[1]*y2[3]+x2[2]*y2[2])+(y2[4]*z2[0]+y2[0]*z2[4]+y2[3]*z2[1]+y2[1]*z2[3]+y2[2]*z2[2])+(z2[4]*x2[0]+z2[0]*x2[4]+z2[3]*x2[1]+z2[1]*x2[3]+z2[2]*x2[2]));
//                    toctic[3]+=2.*((x2[3]*y2[0]+x2[0]*y2[3]+x2[2]*y2[1]+x2[1]*y2[2])+(y2[3]*z2[0]+y2[0]*z2[3]+y2[2]*z2[1]+y2[1]*z2[2])+(z2[3]*x2[0]+z2[0]*x2[3]+z2[2]*x2[1]+z2[1]*x2[2]));
//                    toctic[2]+=2.*((x2[2]*y2[0]+x2[0]*y2[2]+x2[1]*y2[1])+(y2[2]*z2[0]+y2[0]*z2[2]+y2[1]*z2[1])+(z2[2]*x2[0]+z2[0]*x2[2]+z2[1]*x2[1]));
//                    toctic[1]+=2.*((x2[1]*y2[0]+x2[0]*y2[1])+(y2[1]*z2[0]+y2[0]*z2[1])+(z2[1]*x2[0]+z2[0]*x2[1]));
//                    toctic[0]+=2.*((x2[0]*y2[0])+(y2[0]*z2[0])+(z2[0]*x2[0]));
//
//                    toctic[4]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[4]+y2[4])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[4];
//                    toctic[3]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[3]+y2[3])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[3];
//                    toctic[2]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[2]+y2[2])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[2];
//                    toctic[1]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[1]+y2[1])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[1];
//                    toctic[0]+=-2.*(pow(cpocket_r,2.)+pow(ball_radius,2.))*(x2[0]+y2[0])+2.*(pow(cpocket_r,2.)-pow(ball_radius,2.))*z2[0];
//
//                    toctic[0]+=pow(pow(cpocket_r,2.)-pow(ball_radius,2.),2.);
//
//                    for (int k=0;k<8;k++)
//                    {
//                        monoctic[k]=toctic[k]/toctic[8];
//                    }
//                    octic=qsolve_octic(monoctic);
//                    //now to check the roots.
//                    for (int k=0;k<8;k++)
//                    {
//                        if (octic[k]==octic[k] && octic[k]>=0. && octic[k]<t)
//                        {
//                            t0=octic[k];
//                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
//                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
//                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
//                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
//                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
//                            out2[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
//                            a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                            phi=acos(out[2]/a1);
//                            nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                            vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                            vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                            normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                            if (a1<ball_radius && normal>0.)
//                            {
//                                t=t0;
//                                collision=true;
//                            }
//                        }
//                        else
//                        {
//                            c=0;
//                            while (true)
//                            {
//                                c+=1;
//                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0-c*epsilon>=0.)
//                                    {
//                                        t=t0-c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
//                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
//                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
//                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
//                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
//                                out2[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
//                                a1=sqrt(pow(out[2],2.)+pow(cpocket_r-sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.)),2.));
//                                phi=acos(out[2]/a1);
//                                nspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                parspeed=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*sin(atan2(cpockets[j][0]-out[0],cpockets[j][1]-out[1])-atan2(out2[0],out2[1]));
//                                vtheta=parspeed/sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
//                                vphi=sqrt(pow(nspeed,2.)+pow(out2[2],2.))/a1;
//                                normal=gravity*cos(phi)-ball_radius*pow(vphi,2.)+(cpocket_r-ball_radius*sin(phi))*sin(phi)*pow(vtheta,2.);
//                                if (a1<ball_radius && normal>0.)
//                                {
//                                    if (t0+c*epsilon<t)
//                                    {
//                                        t=t0+c*epsilon;
//                                        collision=true;
//                                        break;
//                                    }
//                                }
//                                if (c>1000)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
            }
            else
            {
                for (int j=0;j<2;j++)
                {
                    //middle pockets.
                    t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-mpockets[j][0],2.)+pow(balls[i]._ay[0]-mpockets[j][1],2.)-pow(mpocket_r,2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                            a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                                a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                                a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
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
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]));
                    t0=pow(balls[i]._ax[0]-cpockets[j][0],2.)+pow(balls[i]._ay[0]-cpockets[j][1],2.)-pow(cpocket_r,2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    for (int k=0;k<4;k++)
                    {
                        if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                        {
                            t0=quartic[k];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                            a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                                a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
                                if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                                a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
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
                    t4=pow(balls[i]._ax[2],2.)+pow(balls[i]._ay[2],2.);
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-x)+balls[i]._ay[2]*(balls[i]._ay[0]-y))+pow(balls[i]._ax[1],2.)+pow(balls[i]._ay[1],2.);
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-x)+balls[i]._ay[1]*(balls[i]._ay[0]-y));
                    t0=pow(balls[i]._ax[0]-x,2.)+pow(balls[i]._ay[0]-y,2.)-pow(ctemp[k][2],2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    //check the validity of quartic solutions.
                    for (int l=0;l<4;l++)
                    {
                        if (quartic[l]==quartic[l] && quartic[l]>=0. && quartic[l]<t)
                        {
                            //verify the time of collision.
                            t0=quartic[l];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                            a1=atan2(out[0]-x,out[1]-y);
                            std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                            a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                            a3=ctemp[k][3]+cush[j]._angle;
                            a3=a3-2.*pi*floor(a3/(2.*pi));
                            b3=ctemp[k][4]+cush[j]._angle;
                            b3=b3-2.*pi*floor(b3/(2.*pi));
                            if (a3>pi) {a3=a3-2.*pi;}
                            if (a3<-pi) {a3=a3+2.*pi;}
                            if (b3>pi) {b3=b3-2.*pi;}
                            if (b3<-pi) {b3=b3+2.*pi;}
                            if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
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
                                    out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
                                    {
                                        if (t0-c*epsilon>=0.)
                                        {
                                            t=t0-c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
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

            //check middle straight bits.
            for (int j=0;j<2;j++)
            {
                quadratic=qsolve_quadratic(balls[i]._ax[2],balls[i]._ax[1],balls[i]._ax[0]-mpline[j]);
                for (int k=0;k<2;k++)
                {
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<mpline[0] || out[0]>mpline[1]) && (fabs(blue_y-out[1])>table_width+0.438 && fabs(blue_y-out[1])<table_width+pround*cos(asin(1.625/pround))) && a2>0.)
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

                t4=pow(a1,2.)+pow(b1,2.);
                t3=2.*(a1*a2+b1*b2);
                t2=2.*(a1*a3+b1*b3)+pow(a2,2.)+pow(b2,2.);
                t1=2.*(a2*a3+b2*b3);
                t0=pow(a3,2.)+pow(b3,2.)-pow(ball_radius,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                        b1=cpline[j][0]*a1+cpline[j][1];
                        a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                        a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                        if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                a1=cpline[j][0]*0.5*(out[1]+cpline[j][0]*out[0]-cpline[j][1]);
                                b1=cpline[j][0]*a1+cpline[j][1];
                                a2=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(a1-out[0],b1-out[1])-atan2(out2[0],out2[1]));
                                a3=sqrt(pow(out[0]-a1,2.)+pow(out[1]-b1,2.));
                                if (a1>cpline[j][2] && a1<cpline[j][3] && a2>0. && a3<ball_radius)
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
                a1=sqrt(pow(balls[i]._x-mpockets[j][0],2.)+pow(balls[i]._y-mpockets[j][1],2.));
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
                            a2=sqrt(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                            a3=sqrt(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1],2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]+0.156-0.438-balls[i]._rimpos[k][1]));
                        }
                        else
                        {
                            a2=sqrt(pow(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]-3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                            a3=sqrt(pow(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],2.)+pow(mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1],2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(mpockets[j][0]+3.719-balls[i]._rimpos[k][0],mpockets[j][1]-0.156+0.438-balls[i]._rimpos[k][1]));
                        }
                        if (((a2<ball_radius+2.094 && b2>0.) || (a3<ball_radius+2.094 && b3>0.)))
                        {
                            t=balls[i]._times[k];
                            collision=true;
                            //std::cout << "Hit" << std::endl;
                            break;
                        }
                    }
                    c=j+1;
                    break;
                }
            }
            //std::cout << "Checked middle pockets" << std::endl;
            if (c==0)
            {
                for (int j=0;j<4;j++)
                {
                    a1=sqrt(pow(balls[i]._x-cpockets[j][0],2.)+pow(balls[i]._y-cpockets[j][1],2.));
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
                                a2=sqrt(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==1)
                            {
                                a2=sqrt(pow(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==2)
                            {
                                a2=sqrt(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]-2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]+6.15-balls[i]._rimpos[k][1]));
                            }
                            else if (j==3)
                            {
                                a2=sqrt(pow(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]+2.445-balls[i]._rimpos[k][1],2.));
                                b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]-6.15-balls[i]._rimpos[k][0],cpockets[j][1]+2.445-balls[i]._rimpos[k][1]));
                                a3=sqrt(pow(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],2.)+pow(cpockets[j][1]-6.15-balls[i]._rimpos[k][1],2.));
                                b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(cpockets[j][0]+2.445-balls[i]._rimpos[k][0],cpockets[j][1]-6.15-balls[i]._rimpos[k][1]));
                            }
                            if (((a2<ball_radius+4.5 && b2>0.) || (a3<ball_radius+4.5 && b3>0.)))
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                            //check straight bits.
                            t0=0.5*cpline[2*j][0]*(balls[i]._rimpos[k][1]+cpline[2*j][0]*balls[i]._rimpos[k][0]-cpline[2*j][1]);
                            t1=cpline[2*j][0]*t0+cpline[2*j][1];
                            a2=sqrt(pow(balls[i]._rimpos[k][0]-t0,2.)+pow(balls[i]._rimpos[k][1]-t1,2.));
                            b2=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t0-balls[i]._rimpos[k][0],t1-balls[i]._rimpos[k][1]));

                            t2=0.5*cpline[2*j+1][0]*(balls[i]._rimpos[k][1]+cpline[2*j+1][0]*balls[i]._rimpos[k][0]-cpline[2*j+1][1]);
                            t3=cpline[2*j+1][0]*t2+cpline[2*j+1][1];
                            a3=sqrt(pow(balls[i]._rimpos[k][0]-t2,2.)+pow(balls[i]._rimpos[k][1]-t3,2.));
                            b3=sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-atan2(t2-balls[i]._rimpos[k][0],t3-balls[i]._rimpos[k][1]));
                            if (((a2<ball_radius && t0>cpline[2*j][2] && t0<cpline[2*j][3] && b2>0.) || (a3<ball_radius && t2>cpline[2*j+1][2] && t2<cpline[2*j+1][3] && b3>0.)))
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
                        a1=sqrt(pow(balls[j]._x-mpockets[k][0],2.)+pow(balls[j]._y-mpockets[k][1],2.));
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
                                    a1=sqrt(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2.)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2.));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2.*ball_radius && a2>0.)
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
                        a1=sqrt(pow(balls[j]._x-cpockets[k][0],2.)+pow(balls[j]._y-cpockets[k][1],2.));
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
                                    a1=sqrt(pow(balls[i]._rimpos[m][0]-balls[j]._rimpos[m][0],2.)+pow(balls[i]._rimpos[m][1]-balls[j]._rimpos[m][1],2.));
                                    a2=dot_product(balls[i]._rimvel[m],subtract_vectors(balls[j]._rimpos[m],balls[i]._rimpos[m]))+dot_product(balls[j]._rimvel[m],subtract_vectors(balls[i]._rimpos[m],balls[j]._rimpos[m]));
                                    if (a1<2.*ball_radius && a2>0.)
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
                    xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                    xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                    xmin2-=1.01*ball_radius;
                    xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                    xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                    xmax2+=1.01*ball_radius;
                    ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                    ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                    ymin2-=1.01*ball_radius;
                    ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                    ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                    ymax2+=1.01*ball_radius;
                    //check bounds.
                    if (c<3)
                    {
                        //middle pocket.
                        if (xmin2>blue_x+1.96 || xmax2<blue_x-1.96)
                        {
                            continue;
                        }
                        if ((c==1 && ymax2<mpockets[0][1]-mpocket_r) || (c==2 && ymin2>mpockets[1][1]+mpocket_r))
                        {
                            continue;
                        }
                    }
                    else
                    {
                        //corner pocket.
                        if (c<5 && xmin2>rail_thickness+5.)
                        {
                            continue;
                        }
                        else if (c>4 && xmax2<rail_thickness+2.*cush_thickness+table_length-5.)
                        {
                            continue;
                        }
                        if ((c==3 || c==5) && ymin2>rail_thickness+5.)
                        {
                            continue;
                        }
                        else if ((c==4 || c==6) && ymax2<rail_thickness+2.*cush_thickness+table_width-5.)
                        {
                            continue;
                        }
                    }

                    //potential collision.
                    for (int k=0;k<balls[i]._rimpos.size();k++)
                    {
                        if (balls[i]._times[k]>t)
                        {
                            break;
                        }
                        //check for distance.
                        //assume a 2d collision to approximate.
                        out[0]=balls[j]._ax[2]*balls[i]._times[k]*balls[i]._times[k]+balls[j]._ax[1]*balls[i]._times[k]+balls[j]._ax[0];
                        out[1]=balls[j]._ay[2]*balls[i]._times[k]*balls[i]._times[k]+balls[j]._ay[1]*balls[i]._times[k]+balls[j]._ay[0];
                        if (sqrt(pow(balls[i]._rimpos[k][0]-out[0],2.)+pow(balls[i]._rimpos[k][1]-out[1],2.))-2.*ball_radius<epsilon)
                        {
                            //distance correct.
                            out2[0]=balls[j]._avx[1]*balls[i]._times[k]+balls[j]._avx[0];
                            out2[1]=balls[j]._avy[1]*balls[i]._times[k]+balls[j]._avy[0];
                            b1=atan2(balls[i]._rimpos[k][0]-out[0],balls[i]._rimpos[k][1]-out[1]);
                            a1=sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b1)-sqrt(pow(balls[i]._rimvel[k][0],2.)+pow(balls[i]._rimvel[k][1],2.))*cos(atan2(balls[i]._rimvel[k][0],balls[i]._rimvel[k][1])-b1);
                            if (a1>0.)
                            {
                                t=balls[i]._times[k];
                                collision=true;
                                break;
                            }
                        }
                    }
                }
            }
        }

        //add positions to list for each timestep.
        start=ceil(totaltime/0.01);
        for (int i=0;i<int(floor((totaltime+t)/0.01)-start)+1;i++)
        {
            t0=(start+i)*0.01-totaltime;
            for (int j=0;j<22;j++)
            {
                if (balls[j]._potted)
                {
                    temp[3*j]=-100.;
                    temp[3*j+1]=-100.;
                    temp[3*j+2]=-100.;
                    continue;
                }
                if (!balls[j]._onrim)
                {
                    temp[3*j]=balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0];
                    temp[3*j+1]=balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0];
                    temp[3*j+2]=balls[j]._az[2]*t0*t0+balls[j]._az[1]*t0+balls[j]._az[0];
                }
                else
                {
                    if (int(floor(1000.*t0))==balls[i]._rimpos.size()-1)
                    {
                        out=balls[i]._rimpos[balls[i]._rimpos.size()-1];
                        temp[3*j]=out[0];
                        temp[3*j+1]=out[1];
                        temp[3*j+2]=out[2];
                    }
                    else
                    {
                        //interpolate.
                        a1=balls[j]._times[int(floor(1000.*t0)+1)]-balls[j]._times[int(floor(1000.*t0))];
                        out=balls[j]._rimpos[int(floor(1000.*t0))];
                        out2=balls[j]._rimpos[int(floor(1000.*t0)+1)];
                        temp[3*j]=out[0]+(out2[0]-out[0])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
                        temp[3*j+1]=out[1]+(out2[1]-out[1])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
                        temp[3*j+2]=out[2]+(out2[2]-out[2])*(t0-balls[j]._times[int(floor(1000.*t0))])/a1;
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
                balls[i]._x=balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0];
                balls[i]._y=balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0];
                balls[i]._z=balls[i]._az[2]*t*t+balls[i]._az[1]*t+balls[i]._az[0];
                if (balls[i]._z<-ball_radius)
                {
                    balls[i]._potted=true;
                    balls[i]._vx=0.;
                    balls[i]._vy=0.;
                    balls[i]._vz=0.;
                    balls[i]._xspin=0.;
                    balls[i]._yspin=0.;
                    balls[i]._rspin=0.;
                    pos[pos.size()-1][3*i]=-100.;
                    pos[pos.size()-1][3*i+1]=-100.;
                    pos[pos.size()-1][3*i+2]=-100.;
                    ball_potted_order.push_back(balls[i]._order);
                }
                if (!balls[i]._potted)
                {
                    //velocity.
                    balls[i]._vx=balls[i]._avx[1]*t+balls[i]._avx[0];
                    balls[i]._vy=balls[i]._avy[1]*t+balls[i]._avy[0];
                    balls[i]._vz=balls[i]._avz[1]*t+balls[i]._avz[0];
                    if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                    if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                    if (fabs(balls[i]._vz)<epsilon) {balls[i]._vz=0.;}
                    //spin.
                    balls[i]._xspin=balls[i]._awx[1]*t+balls[i]._awx[0];
                    balls[i]._yspin=balls[i]._awy[1]*t+balls[i]._awy[0];
                    balls[i]._rspin=balls[i]._awz[1]*t+balls[i]._awz[0];
                    if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0.;}
                    if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0.;}
                    if (fabs(balls[i]._rspin)<epsilon) {balls[i]._rspin=0.;}
                }
            }
            else
            {
                //on the rim.
                if (int(floor(1000.*t))==balls[i]._rimpos.size()-1)
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
                    a1=balls[i]._times[int(floor(1000.*t)+1)]-balls[i]._times[int(floor(1000.*t))];
                    out=balls[i]._rimpos[int(floor(1000.*t))];
                    out2=balls[i]._rimpos[int(floor(1000.*t)+1)];
                    balls[i]._x=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._y=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._z=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    out=balls[i]._rimvel[int(floor(1000.*t))];
                    out2=balls[i]._rimvel[int(floor(1000.*t)+1)];
                    balls[i]._vx=out[0]+(out2[0]-out[0])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._vy=out[1]+(out2[1]-out[1])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                    balls[i]._vz=out[2]+(out2[2]-out[2])*(t-balls[i]._times[int(floor(1000.*t))])/a1;
                }
                balls[i]._xspin=balls[i]._vx/ball_radius;
                balls[i]._yspin=balls[i]._vy/ball_radius;
                //work out how rspin changes.
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                if (fabs(balls[i]._vz)<epsilon) {balls[i]._vz=0.;}
                if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0.;}
                if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0.;}
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

//            if (i==5)
//            {
//                std::cout << "x:" << balls[i]._x << std::endl;
//                std::cout << "y:" << balls[i]._y << std::endl;
//                std::cout << "z:" << balls[i]._z << std::endl;
//                std::cout << "vx:" << balls[i]._vx << std::endl;
//                std::cout << "vy:" << balls[i]._vy << std::endl;
//                std::cout << "vz:" << balls[i]._vz << std::endl;
//                std::cout << "xspin:" << double(balls[i]._xspin) << std::endl;
//                std::cout << "yspin:" << double(balls[i]._yspin) << std::endl;
//                std::cout << "rspin:" << double(balls[i]._rspin) << std::endl;
//            }
        }

        //perform collision check.
        if (collision)
        {
            do
            {
                dv=collisions(balls,cush);
                for (int i=0;i<22;i++)
                {
                    if (balls[i]._potted)
                    {
                        continue;
                    }
                    balls[i]._vx+=dv(2*i);
                    balls[i]._vy+=dv(2*i+1);
                    if (fabs(dv(2*i))>=epsilon && fabs(dv(2*i+1))>=epsilon && balls[i]._inflight)
                    {
                        a1=999.;
                        a2=999.;
                        for (int k=0;k<4;k++)
                        {
                            a1=fmin(a1,sqrt(pow(balls[i]._x-cpockets[k][0],2.)+pow(balls[i]._y-cpockets[k][1],2.)));
                        }
                        for (int k=0;k<2;k++)
                        {
                            a2=fmin(a2,sqrt(pow(balls[i]._x-mpockets[k][0],2.)+pow(balls[i]._y-mpockets[k][1],2.)));
                        }
                        if (a1<k1 || (a2<mpocket_r && fabs(blue_y-balls[i]._y)>0.5*table_width+0.438-ball_radius-0.01))
                        {
                            balls[i]._vz=-2*sqrt(pow(balls[i]._vx,2.)+pow(balls[i]._vy,2.));
                        }
                    }
                }
                if (!hitdone)
                {
                    if (fabs(dv(0))>=epsilon || fabs(dv(1))>=epsilon)
                    {
                        for (int j=1;j<22;j++)
                        {
                            if (fabs(dv(2*j))>=epsilon || fabs(dv(2*j+1))>=epsilon)
                            {
                                ball_hit_order.push_back(balls[j]._order);
                            }
                        }
                        hitdone=true;
                    }
                }
            } while (dv!=zero);
            for (int i=0;i<22;i++)
            {
                if (balls[i]._potted)
                {
                    continue;
                }
                if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
                if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
                for (int j=i+1;j<22;j++)
                {
                    t0=sqrt(pow(balls[i]._x-balls[j]._x,2.)+pow(balls[i]._y-balls[j]._y,2.));
                    if (t0-2.*ball_radius<=0.)
                    {
                        t1=atan2(balls[j]._x-balls[i]._x,balls[j]._y-balls[i]._y);
                        t2=(2.*ball_radius-t0)+DOUBLE_EPSILON;
                        balls[i]._x-=t2*sin(t1);
                        balls[i]._y-=t2*cos(t1);
                        balls[j]._x+=t2*sin(t1);
                        balls[j]._y+=t2*cos(t1);
                    }
                }
            }
        }
    }
    return pos;
}

void Server::handleIncomingConnections()
{
    bool alreadyconnected=false;
    for (int i=0;i<2;i++)
    {
        if (players[i].getRemoteAddress()==sf::IpAddress::None)
        {
            //fresh socket ready for connecting to.
            if (listener.accept(players[i])==sf::TcpListener::Done)
            {
                //connected to new socket.
                alreadyconnected=true;
                break;
            }
        }
    }
    if (!alreadyconnected)
    {
        for (int i=0;i<4;i++)
        {
            if (spectators[i].getRemoteAddress()==sf::IpAddress::None)
            {
                //fresh socket ready for connecting to.
                if (listener.accept(spectators[i])==sf::TcpListener::Done)
                {
                    //connected to new socket.
                    break;
                }
            }
        }
    }
}

void Server::executionThread()
{
    sf::Uint16 packetId=0;
    double a=0.;
    double b=0.;
    double c=0.;
    double d=0.;
    double e=0.;
    double f=0.;
    double g=0.;
    std::vector<double> flatresult;
    while (running)
    {
        handleIncomingConnections();

        //check for incoming packets.
        //only receive input from the active player.
        if (players[player_turn].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            if (players[player_turn].receive(packet)==sf::Socket::Done)
            {
                //received an input packet.
                packet >> packetId;
                if (packetId==0)
                {
                    //disconnect this user.
                }
                else if (packetId==1)
                {
                    //trajectory display to other clients.
                    packet >> a >> b >> c >> d >> e >> f >> g;
                    packet.clear();
                    packet << sf::Uint16(0) << a << b << c << d << e << f << g;

                    serverballs[0]._x=f;
                    serverballs[0]._y=g;

                    if (players[!player_turn].getRemoteAddress()!=sf::IpAddress::None)
                    {
                        players[!player_turn].send(packet);
                    }

                    for (int i=0;i<4;i++)
                    {
                        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
                        {
                            spectators[i].send(packet);
                        }
                    }
                }
                else if (packetId==2)
                {
                    //simulate and return the results to everybody.
                    packet >> a >> b >> c >> d >> e;
                    serverballs[0]._vx=a;
                    serverballs[0]._vy=b;
                    serverballs[0]._xspin=c;
                    serverballs[0]._yspin=d;
                    serverballs[0]._rspin=e;
                    result.clear();
                    result=simulate(serverballs,servercushions);

                    respot();

                    placing_white=false;
                    isfoul=false;
                    ismiss=false;
                    ispush=false;
                    foulscore=0;
                    //check if white potted.
                    if (serverballs[0]._potted)
                    {
                        placing_white=true;
                        isfoul=true;
                        foulscore=std::max(foulscore,4);
                        serverballs[0]._potted=false;
                        serverballs[0]._x=cueball_break_x;
                        serverballs[0]._y=cueball_break_y;
                        serverballs[0]._z=ball_radius;
                        serverballs[0]._vx=0.;
                        serverballs[0]._vy=0.;
                        serverballs[0]._vz=0.;
                        serverballs[0]._xspin=0.;
                        serverballs[0]._yspin=0.;
                        serverballs[0]._rspin=0.;
                    }

                    //determine whether or not the shot was legal.
                    if (ball_hit_order.size()==0) {isfoul=true; ismiss=true; foulscore=std::max(foulscore,4);}

                    bool isredon2=false;
                    bool isfreeball2=false;

                    if (isredon)
                    {
                        //pass
                        bool redhit=false;
                        for (int i=0;i<ball_hit_order.size();i++)
                        {
                            if (ball_hit_order[i]<8) {isfoul=true; foulscore=std::max(foulscore,std::max(ball_hit_order[i],4));}
                            if (ball_hit_order[i]>7) {redhit=true;}
                        }
                        if (!redhit) {isfoul=true; ismiss=true; foulscore=std::max(foulscore,4);}
                        for (int i=0;i<ball_potted_order.size();i++)
                        {
                            if (ball_potted_order[i]<8) {isfoul=true; foulscore=std::max(foulscore,std::max(ball_potted_order[i],4));}
                        }
                        if (!isfoul)
                        {
                            //add points.
                            scores[player_turn]+=ball_potted_order.size();
                            current_break+=ball_potted_order.size();
                            if (current_break>highbreak[player_turn]) {highbreak[player_turn]=current_break;}

                            if (ball_potted_order.size()==0)
                            {
                                //no balls potted.
                                current_break=0;
                                player_turn=!player_turn;
                                isredon2=true;
                                isfreeball2=false;
                            }
                            else
                            {
                                isredon2=false;
                                isfreeball2=false;
                            }
                        }
                        else
                        {
                            scores[!player_turn]+=foulscore;
                            current_break=0;
                            player_turn=!player_turn;

                            isfreeball2=is_snookered();
                            if (isfreeball2) {isredon2=false;}
                            else {isredon2=true;}
                        }
                    }
                    else if (!isredon && !isfreeball)
                    {
                        //pass
                        bool colhit=false;
                        for (int i=0;i<ball_hit_order.size();i++)
                        {
                            if (ball_hit_order[i]>7) {isfoul=true; foulscore=std::max(foulscore,4);}
                            if (ball_hit_order[i]<8)
                            {
                                if (ball_hit_order[i]!=nom_colour_order)
                                {
                                    isfoul=true;
                                    foulscore=std::max(foulscore,std::max(ball_hit_order[i],4));
                                }
                                else {colhit=true;}
                            }
                        }
                        if (!colhit) {isfoul=true; ismiss=true; foulscore=std::max(foulscore,4);}
                        bool colpot=false;
                        for (int i=0;i<ball_potted_order.size();i++)
                        {
                            if (ball_potted_order[i]>7) {isfoul=true; foulscore=std::max(foulscore,4);}
                            if (ball_potted_order[i]<8)
                            {
                                if (ball_potted_order[i]!=nom_colour_order)
                                {
                                    isfoul=true;
                                    foulscore=std::max(foulscore,std::max(ball_potted_order[i],4));
                                }
                                else {colpot=true;}
                            }
                        }
                        if (!isfoul)
                        {
                            //add points.
                            if (colpot)
                            {
                                scores[player_turn]+=nom_colour_order;
                                current_break+=nom_colour_order;
                                if (current_break>highbreak[player_turn]) {highbreak[player_turn]=current_break;}
                                isredon2=true;
                                isfreeball2=false;
                            }
                            else
                            {
                                current_break=0;
                                player_turn=!player_turn;
                                isredon2=true;
                                isfreeball2=false;
                            }
                        }
                        else
                        {
                            scores[!player_turn]+=foulscore;
                            current_break=0;
                            player_turn=!player_turn;

                            isfreeball2=is_snookered();
                            if (isfreeball2) {isredon2=false;}
                            else {isredon2=true;}
                        }
                    }
                    else if (!isredon && isfreeball)
                    {
                        //pass
                    }

                    //send the results in a packet.
                    packet.clear();
                    flatresult.clear();
                    for (int i=0;i<result.size();i++)
                    {
                        for (int j=0;j<66;j++)
                        {
                            flatresult.push_back(result[i][j]);
                        }
                    }
                    for (int i=0;i<22;i++)
                    {
                        if (!serverballs[i]._potted)
                        {
                            flatresult.push_back(serverballs[i]._x);
                            flatresult.push_back(serverballs[i]._y);
                            flatresult.push_back(serverballs[i]._z);
                        }
                        else
                        {
                            flatresult.push_back(-100.);
                            flatresult.push_back(-100.);
                            flatresult.push_back(-100.);
                        }
                    }
                    packet << sf::Uint16(1) << sf::Uint32(result.size()+1);
                    for (int i=0;i<flatresult.size();i++)
                    {
                        packet << flatresult[i];
                    }

                    for (int i=0;i<2;i++)
                    {
                        if (players[i].getRemoteAddress()!=sf::IpAddress::None)
                        {
                            players[i].send(packet);
                        }
                    }
                    for (int i=0;i<4;i++)
                    {
                        if (spectators[i].getRemoteAddress()!=sf::IpAddress::None)
                        {
                            spectators[i].send(packet);
                        }
                    }
                    turnpacket();
                }
                else if (packetId==3)
                {
                    //place again after foul.
                }
                else if (packetId==4)
                {
                    //concede frame.
                    frames[!player_turn]+=1;
                    resetframe();
                    turnpacket();
                }
                else if (packetId==5)
                {
                    //concede match.
                    frames[!player_turn]=(framesbestof+1)/2;
                    gameover=true;
                    turnpacket();
                }
            }
        }

        if (players[!player_turn].getRemoteAddress()!=sf::IpAddress::None)
        {
            packet.clear();
            if (players[!player_turn].receive(packet)==sf::Socket::Done)
            {
                packet >> packetId;
                if (packetId==4)
                {
                    //concede frame.
                    frames[player_turn]+=1;
                    resetframe();
                    turnpacket();
                }
                else if (packetId==5)
                {
                    //concede match.
                    frames[player_turn]=(framesbestof+1)/2;
                    gameover=true;
                    turnpacket();
                }
            }
        }

        sf::sleep(sf::milliseconds(5));
    }
}

Eigen::Matrix<double,46,1> collisions2(Ball b[22],Cushion cush[6])
{
    //broad phase collision check.
    Eigen::MatrixXd dv=Eigen::MatrixXd::Zero(46,1);

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
            if (b[i]._vx!=0. || b[i]._vy!=0. || b[i]._vz!=0.)
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
                                    dist=sqrt(pow(dx,2.)+pow(dy,2.));
                                    relspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(dx,dy)-atan2(b[i]._vx,b[i]._vy))-sqrt(pow(b[thing-1]._vx,2.)+pow(b[thing-1]._vy,2.))*cos(atan2(dx,dy)-atan2(b[thing-1]._vx,b[thing-1]._vy));

                                    if (dist-2.*ball_radius<epsilon && relspeed>0.)
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
                                        col+=1;
                                    }
                                }
                            }
                        }
                    }
                }
                //check cushion collisions.
                if (fabs(blue_x-b[i]._x)>table_width-ball_radius-2.*epsilon || fabs(blue_y-b[i]._y)>0.5*table_width-ball_radius-2.*epsilon)
                {
                    double min_dist=999.;
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
                    relspeed=-sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(normal_angle-atan2(b[i]._vx,b[i]._vy));
                    if (min_dist-ball_radius<epsilon && relspeed>0.)
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
                        if (relspin>0.)
                        {
                            //topspin.
                            dw=-5.*relspeed*(muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin<0.)
                        {
                            //backspin.
                            dw=-5.*relspeed*(-muk*ball_radius*cos(cushion_alpha)+cushion_diff)/pow(ball_radius,2.);
                        }
                        if (relspin==0.)
                        {
                            dw=-5.*relspeed*cushion_diff/pow(ball_radius,2.);
                        }
                        dw=fmax(dw,-relspin-relspeed/ball_radius);
                        b[i]._xspin+=dw*sin(normal_angle+pi);
                        b[i]._yspin+=dw*cos(normal_angle+pi);
                        //adjust sidespin off cushion.
                        double parspeed=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*sin(atan2(b[i]._vx,b[i]._vy)-normal_angle);
                        double dvpar;

                        if ((parspeed>=0. && b[i]._rspin*ball_radius>parspeed) || (parspeed<0. && b[i]._rspin*ball_radius>parspeed))
                        {
                            dw=-5.*muk*relspeed*cos(cushion_alpha)/ball_radius;
                            dvpar=-0.4*dw*ball_radius;
                            if ((b[i]._rspin+dw)*ball_radius<parspeed+dvpar)
                            {
                                dw=-5.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                                dvpar=2.*(b[i]._rspin*ball_radius-parspeed)/(7.*ball_radius);
                            }
                        }
                        else if ((parspeed>=0. && b[i]._rspin*ball_radius<parspeed) || (parspeed<0. && b[i]._rspin*ball_radius<parspeed))
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
                        std::tie(_speed,_angle)=add_vectors(sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.)),atan2(b[i]._vx,b[i]._vy),dvpar,normal_angle+0.5*pi);
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
                    std::tie(v,an)=add_vectors(-sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(b[i]._vx,b[i]._vy)-atan2(Z(2*i,0),Z(2*i+1,0))),atan2(Z(2*i,0),Z(2*i+1,0)),sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(b[i]._vx,b[i]._vy)-(atan2(Z(2*i,0),Z(2*i+1,0))+0.5*pi)),atan2(Z(2*i,0),Z(2*i+1,0))+0.5*pi);
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
                            double pari=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(b[i]._vx,b[i]._vy)-na);
                            double peri=sqrt(pow(b[i]._vx,2.)+pow(b[i]._vy,2.))*cos(atan2(b[i]._vx,b[i]._vy)-(na+0.5*pi));
                            double parj=sqrt(pow(b[j]._vx,2.)+pow(b[j]._vy,2.))*cos(atan2(b[j]._vx,b[j]._vy)-na);
                            double perj=sqrt(pow(b[j]._vx,2.)+pow(b[j]._vy,2.))*cos(atan2(b[j]._vx,b[j]._vy)-(na+0.5*pi));
                            double speed;
                            double theta;
                            std::tie(speed,theta)=add_vectors(parj,na,peri,na+0.5*pi);
                            dv(2*i)=speed*sin(theta)-b[i]._vx;
                            dv(2*i+1)=speed*cos(theta)-b[i]._vy;
                            std::tie(speed,theta)=add_vectors(pari,na,perj,na+0.5*pi);
                            dv(2*j)=speed*sin(theta)-b[j]._vx;
                            dv(2*j+1)=speed*cos(theta)-b[j]._vy;
                            //adjust spins of balls.
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    //simultaneous.
    if (col>0)
    {
        Eigen::MatrixXd v(46,1);
        Eigen::MatrixXd J(col,1);

        for (int i=0;i<22;i++)
        {
            v(2*i,0)=b[i]._vx;
            v(2*i+1,0)=b[i]._vy;
        }
        v(44,0)=0;
        v(45,0)=0;

        //J=(Z.transpose()*M_*Z).colPivHouseholderQr().solve(-(1.+eball)*Z.transpose()*v);
        J=(Z.transpose()*M_*Z).householderQr().solve(-(1.+eball)*Z.transpose()*v);
        dv=M_*Z*J;
    }
    return dv;
}

std::array<std::vector<double>,3> trajectory(Ball balls[22],Cushion cush[6])
{
    //assume the only moving ball is the white and plot its trajectory.
    static double t_object=1.5;
    static double t_cueball=3.;

    double t=0.;
    double totaltime=0.;
    double start=0.;
    int c=0;

    double tstart_object=0.;
    double tstart_cueball=0.;

    std::array<std::vector<double>,3> pos;
    std::array<double,6> temp;
    std::array<double,3> out;
    std::array<double,3> out2;
    std::array<double,3> out3;
    std::array<double,3> out4;
    std::array<double,4> quartic;
    std::array<double,2> quadratic;
    std::array<std::array<double,5>,6> ctemp;
    Eigen::Matrix<double,46,1> dv;
    Eigen::Matrix<double,46,1> zero;
    std::vector<int> indices;
    indices.push_back(0);

    for (int i=0;i<46;i++)
    {
        dv(i,0)=0.;
        zero(i,0)=0.;
    }

    double a1=0.;
    double b1=0.;
    double c1=0.;
    double a2=0.;
    double b2=0.;
    double c2=0.;
    double a3=0.;
    double b3=0.;
    double c3=0.;
    double t0=0.;
    double t1=0.;
    double t2=0.;
    double t3=0.;
    double t4=0.;
    double x=0.;
    double y=0.;

    double xmin=0.;
    double xmax=0.;
    double xmin2=0.;
    double xmax2=0.;
    double ymin=0.;
    double ymax=0.;
    double ymin2=0.;
    double ymax2=0.;

    int xpos=0;
    int ypos=0;
    int gxmin=0;
    int gxmax=0;
    int gymin=0;
    int gymax=0;

    bool collision=false;
    int i=0;
    int j=0;
    int object=0;

    bool cueball_done1=false;
    bool cueball_done2=false;
    bool object_done=false;

    //save initial positions.
    double originalpos[2][2]={};

    originalpos[0][0]=balls[0]._x;
    originalpos[0][1]=balls[0]._y;
    originalpos[1][0]=balls[0]._x;
    originalpos[1][1]=balls[0]._y;

    for (int i=0;i<22;i++)
    {
        if (balls[i]._x>0)
        {
            balls[i].update_equation();
        }
    }

    while (true)
    {
        collision=false;
        t=999.;

        //std::cout << "Time: " << totaltime << std::endl;
//        if (totaltime>1.97)
//        {
//            break;
//        }

        for (int abc=0;abc<indices.size();abc++)
        {
            i=indices[abc];
            balls[i].update_equation();
            t=fmin(t,balls[i]._t);
        }

        if (balls[0]._onrim || balls[0]._inflight)
        {
            if (!cueball_done1)
            {
                break;
            }
            else
            {
                cueball_done2=true;
                indices.clear();
                if (!object_done)
                {
                    indices.push_back(object);
                }
            }
        }
        if (object!=0 && (balls[object]._onrim || balls[object]._inflight))
        {
            indices.clear();
            if (!cueball_done2)
            {
                indices.push_back(0);
            }
            object_done=true;
        }

        if (t>998.)
        {
            break;
        }

        for (int abc=0;abc<indices.size();abc++)
        {
            i=indices[abc];

            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmin-=1.01*ball_radius;
            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
            xmax+=1.01*ball_radius;
            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymin-=1.01*ball_radius;
            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
            ymax+=1.01*ball_radius;

            for (int j=1;j<22;j++)
            {
                if (balls[j]._x<0 || i==j || (i==0 && j==object))
                {
                    continue;
                }

                //initial check for i.
                xmin2=fmin(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmin2=fmin(xmin2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmin2-=1.01*ball_radius;
                xmax2=fmax(balls[j]._x,balls[j]._x-(balls[j]._ax[1]*balls[j]._ax[1])/(4.*balls[j]._ax[2]));
                xmax2=fmax(xmax2,balls[j]._ax[2]*t*t+balls[j]._ax[1]*t+balls[j]._ax[0]);
                xmax2+=1.01*ball_radius;
                ymin2=fmin(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymin2=fmin(ymin2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymin2-=1.01*ball_radius;
                ymax2=fmax(balls[j]._y,balls[j]._y-(balls[j]._ay[1]*balls[j]._ay[1])/(4.*balls[j]._ay[2]));
                ymax2=fmax(ymax2,balls[j]._ay[2]*t*t+balls[j]._ay[1]*t+balls[j]._ay[0]);
                ymax2+=1.01*ball_radius;

                if (!((xmax>xmin2 || xmax2>xmin) && (ymax>ymin2 || ymax2>ymin)))
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

                t0=pow(c1,2.)+pow(c2,2.)+pow(c3,2.)-4.*pow(ball_radius,2.);
                t1=2.*(b1*c1+b2*c2+b3*c3);
                t2=2.*(a1*c1+a2*c2+a3*c3)+pow(b1,2.)+pow(b2,2.)+pow(b3,2.);
                t3=2.*(a1*b1+a2*b2+a3*b3);
                t4=pow(a1,2.)+pow(a2,2.)+pow(a3,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        //verify the time of actual collision to be certain.
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                        out2[0]=balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0];
                        out2[1]=balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0];
                        out2[2]=balls[j]._az[2]*t0*t0+balls[j]._az[1]*t0+balls[j]._az[0];
                        out3[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        out3[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        out3[2]=balls[i]._avz[1]*t0+balls[i]._avz[0];
                        out4[0]=balls[j]._avx[1]*t0+balls[j]._avx[0];
                        out4[1]=balls[j]._avy[1]*t0+balls[j]._avy[0];
                        out4[2]=balls[j]._avz[1]*t0+balls[j]._avz[0];
                        a1=get_relspeed(out,out3,out2,out4);
                        if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                        {
                            //collision at the specified time.
                            t=t0;
                            collision=true;

                            xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmin-=1.01*ball_radius;
                            xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                            xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                            xmax+=1.01*ball_radius;
                            ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymin-=1.01*ball_radius;
                            ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                            ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                            ymax+=1.01*ball_radius;
                        }
                        else
                        {
                            //adjust the time minutely to ensure a collision.
                            c=0;
                            while (true)
                            {
                                c+=1;
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ax[1]*(t0-c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._ay[1]*(t0-c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[j]._az[1]*(t0-c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0-c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0-c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0-c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0-c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                out2[0]=balls[j]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ax[1]*(t0+c*epsilon)+balls[j]._ax[0];
                                out2[1]=balls[j]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._ay[1]*(t0+c*epsilon)+balls[j]._ay[0];
                                out2[2]=balls[j]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[j]._az[1]*(t0+c*epsilon)+balls[j]._az[0];
                                out3[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                out3[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                out3[2]=balls[i]._avz[1]*(t0+c*epsilon)+balls[i]._avz[0];
                                out4[0]=balls[j]._avx[1]*(t0+c*epsilon)+balls[j]._avx[0];
                                out4[1]=balls[j]._avy[1]*(t0+c*epsilon)+balls[j]._avy[0];
                                out4[2]=balls[j]._avz[1]*(t0+c*epsilon)+balls[j]._avz[0];
                                a1=get_relspeed(out,out3,out2,out4);
                                if (sqrt(pow(out[0]-out2[0],2.)+pow(out[1]-out2[1],2.)+pow(out[2]-out2[2],2.))-2.*ball_radius<epsilon && a1>0.)
                                {
                                    if(t0+c*epsilon<t)
                                    {
                                        t=t0+c*epsilon;
                                        collision=true;

                                        xmin=fmin(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmin=fmin(xmin,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmin-=1.01*ball_radius;
                                        xmax=fmax(balls[i]._x,balls[i]._x-(balls[i]._ax[1]*balls[i]._ax[1])/(4.*balls[i]._ax[2]));
                                        xmax=fmax(xmax,balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0]);
                                        xmax+=1.01*ball_radius;
                                        ymin=fmin(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymin=fmin(ymin,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymin-=1.01*ball_radius;
                                        ymax=fmax(balls[i]._y,balls[i]._y-(balls[i]._ay[1]*balls[i]._ay[1])/(4.*balls[i]._ay[2]));
                                        ymax=fmax(ymax,balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0]);
                                        ymax+=1.01*ball_radius;
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
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                        if (j==0) {a2=-out2[0];}
                        else {a2=out2[0];}
                        if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                if (j==0) {a2=-out2[0];}
                                else {a2=out2[0];}
                                if ((out[0]<xlim[0] || out[0]>xlim[1]) && (out[1]>rail_thickness+5.43 && out[1]<rail_thickness+2.*cush_thickness+table_width-5.43) && a2>0.)
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
                    if (quadratic[k]==quadratic[k] && quadratic[k]>=0. && quadratic[k]<t)
                    {
                        t0=quadratic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                        if (j==0) {a2=-out2[1];}
                        else {a2=out2[1];}
                        if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
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
                                out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
                                {
                                    if (t0-c*epsilon>=0.)
                                    {
                                        t=t0-c*epsilon;
                                        collision=true;
                                        break;
                                    }
                                }
                                out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                if (j==0) {a2=-out2[1];}
                                else {a2=out2[1];}
                                if ((out[1]<ylim[0] || out[1]>ylim[1]) && ((out[0]>rail_thickness+5.43 && out[0]<rail_thickness+cush_thickness+table_width-4.344) || (out[0]>rail_thickness+cush_thickness+table_width+4.344 && out[0]<rail_thickness+2.*cush_thickness+table_length-5.43)) && a2>0.)
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

            for (int j=0;j<2;j++)
            {
                //middle pockets.
                t4=balls[i]._ax[2]*balls[i]._ax[2]+balls[i]._ay[2]*balls[i]._ay[2];
                t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-mpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-mpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-mpockets[j][1]));
                t0=pow(balls[i]._ax[0]-mpockets[j][0],2.)+pow(balls[i]._ay[0]-mpockets[j][1],2.)-pow(mpocket_r,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                        a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
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
                            out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                            a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
                            if (a1<mpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                if (t0-c*epsilon>=0.)
                                {
                                    t=t0-c*epsilon;
                                    break;
                                }
                            }
                            out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-mpockets[j][0],2.)+pow(out[1]-mpockets[j][1],2.));
                            a2=sqrt(pow(mpocket_r-a1,2.)+pow(ball_radius,2.));
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
                t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[2]*(balls[i]._ay[0]-cpockets[j][1]))+balls[i]._ax[1]*balls[i]._ax[1]+balls[i]._ay[1]*balls[i]._ay[1];
                t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-cpockets[j][0])+balls[i]._ay[1]*(balls[i]._ay[0]-cpockets[j][1]));
                t0=pow(balls[i]._ax[0]-cpockets[j][0],2.)+pow(balls[i]._ay[0]-cpockets[j][1],2.)-pow(cpocket_r,2.);

                quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                for (int k=0;k<4;k++)
                {
                    if (quartic[k]==quartic[k] && quartic[k]>=0. && quartic[k]<t)
                    {
                        t0=quartic[k];
                        out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                        out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                        a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                        a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
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
                            out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                            a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
                            if (a1<cpocket_r+epsilon && a2<=ball_radius+epsilon)
                            {
                                if (t0-c*epsilon>=0.)
                                {
                                    t=t0-c*epsilon;
                                    break;
                                }
                            }
                            out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                            a1=sqrt(pow(out[0]-cpockets[j][0],2.)+pow(out[1]-cpockets[j][1],2.));
                            a2=sqrt(pow(cpocket_r-a1,2.)+pow(ball_radius,2.));
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
                    t4=pow(balls[i]._ax[2],2.)+pow(balls[i]._ay[2],2.);
                    t3=2.*(balls[i]._ax[2]*balls[i]._ax[1]+balls[i]._ay[2]*balls[i]._ay[1]);
                    t2=2.*(balls[i]._ax[2]*(balls[i]._ax[0]-x)+balls[i]._ay[2]*(balls[i]._ay[0]-y))+pow(balls[i]._ax[1],2.)+pow(balls[i]._ay[1],2.);
                    t1=2.*(balls[i]._ax[1]*(balls[i]._ax[0]-x)+balls[i]._ay[1]*(balls[i]._ay[0]-y));
                    t0=pow(balls[i]._ax[0]-x,2.)+pow(balls[i]._ay[0]-y,2.)-pow(ctemp[k][2],2.);

                    quartic=qsolve_quartic(t4,t3,t2,t1,t0);

                    //check the validity of quartic solutions.
                    for (int l=0;l<4;l++)
                    {
                        if (quartic[l]==quartic[l] && quartic[l]>=0. && quartic[l]<t)
                        {
                            //verify the time of collision.
                            t0=quartic[l];
                            out[0]=balls[i]._ax[2]*t0*t0+balls[i]._ax[1]*t0+balls[i]._ax[0];
                            out[1]=balls[i]._ay[2]*t0*t0+balls[i]._ay[1]*t0+balls[i]._ay[0];
                            out[2]=balls[i]._az[2]*t0*t0+balls[i]._az[1]*t0+balls[i]._az[0];
                            out2[0]=balls[i]._avx[1]*t0+balls[i]._avx[0];
                            out2[1]=balls[i]._avy[1]*t0+balls[i]._avy[0];
                            a1=atan2(out[0]-x,out[1]-y);
                            std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                            a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                            a3=ctemp[k][3]+cush[j]._angle;
                            a3=a3-2.*pi*floor(a3/(2.*pi));
                            b3=ctemp[k][4]+cush[j]._angle;
                            b3=b3-2.*pi*floor(b3/(2.*pi));
                            if (a3>pi) {a3=a3-2.*pi;}
                            if (a3<-pi) {a3=a3+2.*pi;}
                            if (b3>pi) {b3=b3-2.*pi;}
                            if (b3<-pi) {b3=b3+2.*pi;}

                            if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
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
                                    out[0]=balls[i]._ax[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ax[1]*(t0-c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._ay[1]*(t0-c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0-c*epsilon)*(t0-c*epsilon)+balls[i]._az[1]*(t0-c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0-c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0-c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
                                    {
                                        if (t0-c*epsilon>=0.)
                                        {
                                            t=t0-c*epsilon;
                                            collision=true;
                                            break;
                                        }
                                    }
                                    out[0]=balls[i]._ax[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ax[1]*(t0+c*epsilon)+balls[i]._ax[0];
                                    out[1]=balls[i]._ay[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._ay[1]*(t0+c*epsilon)+balls[i]._ay[0];
                                    out[2]=balls[i]._az[2]*(t0+c*epsilon)*(t0+c*epsilon)+balls[i]._az[1]*(t0+c*epsilon)+balls[i]._az[0];
                                    out2[0]=balls[i]._avx[1]*(t0+c*epsilon)+balls[i]._avx[0];
                                    out2[1]=balls[i]._avy[1]*(t0+c*epsilon)+balls[i]._avy[0];
                                    a1=atan2(out[0]-x,out[1]-y);
                                    std::tie(b1,b2)=cush[j].distance(out[0],out[1],out[2]);
                                    a2=-sqrt(pow(out2[0],2.)+pow(out2[1],2.))*cos(atan2(out2[0],out2[1])-b2);
                                    if (((a3<b3 && a1>a3 && a1<b3) || (a3>b3 && (a1>a3 || a1<b3))) && b1-ball_radius<epsilon && a2>0.)
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
        }

        //add positions to list for each timestep.
        start=ceil(totaltime/0.01);
        for (int i=0;i<int(floor((totaltime+t)/0.01)-start)+1;i++)
        {
            t0=(start+i)*0.01-totaltime;
            if (cueball_done1 && totaltime+t0>tstart_cueball+t_cueball)
            {
                cueball_done2=true;
                indices.clear();
                if (object!=0 && !object_done)
                {
                    indices.push_back(object);
                }
            }
            if (cueball_done1 && totaltime+t0>tstart_object+t_object)
            {
                object_done=true;
                indices.clear();
                if (!cueball_done2)
                {
                    indices.push_back(0);
                }
            }

            if (!cueball_done1)
            {
                pos[0].push_back(balls[0]._ax[2]*t0*t0+balls[0]._ax[1]*t0+balls[0]._ax[0]);
                pos[0].push_back(balls[0]._ay[2]*t0*t0+balls[0]._ay[1]*t0+balls[0]._ay[0]);
            }
            else
            {
                if (!cueball_done2)
                {
                    pos[1].push_back(balls[0]._ax[2]*t0*t0+balls[0]._ax[1]*t0+balls[0]._ax[0]);
                    pos[1].push_back(balls[0]._ay[2]*t0*t0+balls[0]._ay[1]*t0+balls[0]._ay[0]);
                }
                if (!object_done)
                {
                    j=indices[1];
                    pos[2].push_back(balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0]);
                    pos[2].push_back(balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0]);
                }
            }
        }
        t0=t;
        if (cueball_done1 && totaltime+t0>tstart_cueball+t_cueball)
        {
            cueball_done2=true;
            indices.clear();
            if (object!=0 && !object_done)
            {
                indices.push_back(object);
            }
        }
        if (cueball_done1 && totaltime+t0>tstart_object+t_object)
        {
            object_done=true;
            indices.clear();
            if (!cueball_done2)
            {
                indices.push_back(0);
            }
        }

        if (!cueball_done1)
        {
            pos[0].push_back(balls[0]._ax[2]*t0*t0+balls[0]._ax[1]*t0+balls[0]._ax[0]);
            pos[0].push_back(balls[0]._ay[2]*t0*t0+balls[0]._ay[1]*t0+balls[0]._ay[0]);
        }
        else
        {
            if (!cueball_done2)
            {
                pos[1].push_back(balls[0]._ax[2]*t0*t0+balls[0]._ax[1]*t0+balls[0]._ax[0]);
                pos[1].push_back(balls[0]._ay[2]*t0*t0+balls[0]._ay[1]*t0+balls[0]._ay[0]);
            }
            if (!object_done)
            {
                j=indices[1];
                pos[2].push_back(balls[j]._ax[2]*t0*t0+balls[j]._ax[1]*t0+balls[j]._ax[0]);
                pos[2].push_back(balls[j]._ay[2]*t0*t0+balls[j]._ay[1]*t0+balls[j]._ay[0]);
            }
        }

        totaltime+=t;

        //check if done because of time.
        if (object_done && cueball_done2 && cueball_done1)
        {
            break;
        }

        grid={};
        grid_index={};
        for (int abc=0;abc<indices.size();abc++)
        {
            i=indices[abc];
            balls[i]._x=balls[i]._ax[2]*t*t+balls[i]._ax[1]*t+balls[i]._ax[0];
            balls[i]._y=balls[i]._ay[2]*t*t+balls[i]._ay[1]*t+balls[i]._ay[0];
            //velocity.
            balls[i]._vx=balls[i]._avx[1]*t+balls[i]._avx[0];
            balls[i]._vy=balls[i]._avy[1]*t+balls[i]._avy[0];
            if (fabs(balls[i]._vx)<epsilon) {balls[i]._vx=0.;}
            if (fabs(balls[i]._vy)<epsilon) {balls[i]._vy=0.;}
            //spin.
            balls[i]._xspin=balls[i]._awx[1]*t+balls[i]._awx[0];
            balls[i]._yspin=balls[i]._awy[1]*t+balls[i]._awy[0];
            balls[i]._rspin=balls[i]._awz[1]*t+balls[i]._awz[0];
            if (fabs(balls[i]._xspin)<epsilon) {balls[i]._xspin=0.;}
            if (fabs(balls[i]._yspin)<epsilon) {balls[i]._yspin=0.;}
            if (fabs(balls[i]._rspin)<epsilon) {balls[i]._rspin=0.;}
        }
        for (int i=0;i<22;i++)
        {
            if (balls[i]._x<0)
            {
                continue;
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
        }

        if (collision)
        {
            if (!cueball_done1)
            {
                cueball_done1=true;
                tstart_cueball=totaltime;
                tstart_object=totaltime;
                //first collision - do it.
                //afterwards make all balls except for the object ball and cueball still.
                do
                {
                    dv=collisions2(balls,cush);
                    for (int i=0;i<22;i++)
                    {
                        balls[i]._vx+=dv(2*i);
                        balls[i]._vy+=dv(2*i+1);
                    }
                } while (dv!=zero);

                for (int j=1;j<22;j++)
                {
                    a1=sqrt(pow(balls[j]._x-balls[0]._x,2.)+pow(balls[j]._y-balls[0]._y,2.));
                    if (a1-2.*ball_radius<epsilon && object==0)
                    {
                        //move the balls apart.
                        object=j;
                        indices.push_back(object);
                        originalpos[1][0]=balls[object]._x;
                        originalpos[1][1]=balls[object]._y;
                        if (a1-2.*ball_radius<=0.)
                        {
                            t1=atan2(balls[j]._x-balls[0]._x,balls[j]._y-balls[0]._y);
                            t2=(2.*ball_radius-a1)+DOUBLE_EPSILON;
                            balls[0]._x-=t2*sin(t1);
                            balls[0]._y-=t2*cos(t1);
                            balls[j]._x+=t2*sin(t1);
                            balls[j]._y+=t2*cos(t1);
                        }
                    }
                    else
                    {
                        balls[j]._vx=0.;
                        balls[j]._vy=0.;
                    }
                }
                if (object==0)
                {
                    //hits a cushion.
                    break;
                }
            }
            else
            {
                //check which ball is at a new event.
                if (!cueball_done2)
                {
                    for (int j=1;j<22;j++)
                    {
                        //check cueball.
                        a1=sqrt(pow(balls[0]._x-balls[j]._x,2.)+pow(balls[0]._y-balls[j]._y,2.));
                        a2=sqrt(pow(balls[0]._vx,2.)+pow(balls[0]._vy,2.))*cos(atan2(balls[j]._x-balls[0]._x,balls[j]._y-balls[0]._y)-atan2(balls[0]._vx,balls[0]._vy));
                        a2-=sqrt(pow(balls[j]._vx,2.)+pow(balls[j]._vy,2.))*cos(atan2(balls[j]._x-balls[0]._x,balls[j]._y-balls[0]._y)-atan2(balls[j]._vx,balls[j]._vy));
                        if (a1-2.*ball_radius<epsilon && a2>0.)
                        {
                            cueball_done2=true;
                            balls[0]._vx=0.;
                            balls[0]._vy=0.;
                            balls[0]._xspin=0.;
                            balls[0]._yspin=0.;
                            balls[0]._rspin=0.;
                            indices.clear();
                            if (!object_done && object!=0)
                            {
                                indices.push_back(object);
                            }
//                            if (a1-2.*ball_radius<=0.)
//                            {
//                                t1=atan2(balls[j]._x-balls[0]._x,balls[j]._y-balls[0]._y);
//                                t2=(2.*ball_radius-a1)+DOUBLE_EPSILON;
//                                balls[0]._x-=t2*sin(t1);
//                                balls[0]._y-=t2*cos(t1);
//                                balls[j]._x+=t2*sin(t1);
//                                balls[j]._y+=t2*cos(t1);
//                            }
                            break;
                        }
                    }
                    for (int j=0;j<6;j++)
                    {
                        std::tie(a1,a2)=cush[j].distance(balls[0]._x,balls[0]._y,ball_radius);
                        b2=-sqrt(pow(balls[0]._vx,2.)+pow(balls[0]._vy,2.))*cos(atan2(balls[0]._vx,balls[0]._vy)-a2);

                        //std::cout << "Distance to cushion: " << a1 << std::endl;
                        if (a1-ball_radius<epsilon && b2>0.)
                        {
                            //hits a cushion.
                            cueball_done2=true;
                            balls[0]._vx=0.;
                            balls[0]._vy=0.;
                            balls[0]._xspin=0.;
                            balls[0]._yspin=0.;
                            balls[0]._rspin=0.;
                            indices.clear();
                            if (!object_done && object!=0)
                            {
                                indices.push_back(object);
                            }
                            break;
                        }
                    }
                }
                if (!object_done)
                {
                    for (int j=0;j<22;j++)
                    {
                        if (j!=object)
                        {
                            a1=sqrt(pow(balls[object]._x-balls[j]._x,2.)+pow(balls[object]._y-balls[j]._y,2.));
                            a2=sqrt(pow(balls[object]._vx,2.)+pow(balls[object]._vy,2.))*cos(atan2(balls[j]._x-balls[object]._x,balls[j]._y-balls[object]._y)-atan2(balls[object]._vx,balls[object]._vy));
                            a2-=sqrt(pow(balls[j]._vx,2.)+pow(balls[j]._vy,2.))*cos(atan2(balls[j]._x-balls[object]._x,balls[j]._y-balls[object]._y)-atan2(balls[j]._vx,balls[j]._vy));
                            if (a1-2.*ball_radius<epsilon && a2>0.)
                            {
                                object_done=true;
                                balls[object]._vx=0.;
                                balls[object]._vy=0.;
                                balls[object]._xspin=0.;
                                balls[object]._yspin=0.;
                                balls[object]._rspin=0.;
                                indices.clear();
                                if (!cueball_done2)
                                {
                                    indices.push_back(0);
                                }
                                //std::cout << "INDICES SIZE" << indices.size() << std::endl;
//                                if (a1-2.*ball_radius<=0.)
//                                {
//                                    t1=atan2(balls[j]._x-balls[object]._x,balls[j]._y-balls[object]._y);
//                                    t2=(2.*ball_radius-a1)+DOUBLE_EPSILON;
//                                    balls[object]._x-=t2*sin(t1);
//                                    balls[object]._y-=t2*cos(t1);
//                                    balls[j]._x+=t2*sin(t1);
//                                    balls[j]._y+=t2*cos(t1);
//                                }
                                break;
                            }
                        }
                    }
                    for (int j=0;j<6;j++)
                    {
                        std::tie(a1,a2)=cush[j].distance(balls[object]._x,balls[object]._y,ball_radius);
                        b2=-sqrt(pow(balls[object]._vx,2.)+pow(balls[object]._vy,2.))*cos(atan2(balls[object]._vx,balls[object]._vy)-a2);
                        if (a1-ball_radius<epsilon && b2>0.)
                        {
                            //hits a cushion.
                            object_done=true;
                            balls[object]._vx=0.;
                            balls[object]._vy=0.;
                            balls[object]._xspin=0.;
                            balls[object]._yspin=0.;
                            balls[object]._rspin=0.;
                            indices.clear();
                            if (!cueball_done2)
                            {
                                indices.push_back(0);
                            }
                            //std::cout << "INDICES SIZE" << indices.size() << std::endl;
                            break;
                        }
                    }
                }
            }
        }

        if (object_done && cueball_done2 && cueball_done1)
        {
            break;
        }
    }
    //change the positions back.
    balls[0]._x=originalpos[0][0];
    balls[0]._y=originalpos[0][1];
    balls[0]._vx=0.;
    balls[0]._vy=0.;
    balls[0]._xspin=0.;
    balls[0]._yspin=0.;
    balls[0]._rspin=0.;
    balls[object]._x=originalpos[1][0];
    balls[object]._y=originalpos[1][1];
    balls[object]._vx=0.;
    balls[object]._vy=0.;
    balls[object]._xspin=0.;
    balls[object]._yspin=0.;
    balls[object]._rspin=0.;

    return pos;
}

#endif // OBJECTS_H_INCLUDED
