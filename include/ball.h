#ifndef BALL_H_INCLUDED
#define BALL_H_INCLUDED

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
    _shape.setOrigin(dfactor*_r,dfactor*_r);
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

#endif // BALL_H_INCLUDED
