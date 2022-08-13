#ifndef BALL-TRAJECTORY_H_INCLUDED
#define BALL-TRAJECTORY_H_INCLUDED

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

#endif // BALL-TRAJECTORY_H_INCLUDED
