#ifndef CUSHION_H_INCLUDED
#define CUSHION_H_INCLUDED

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

#endif // CUSHION_H_INCLUDED
