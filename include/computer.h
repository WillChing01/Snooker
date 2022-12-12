#ifndef COMPUTER_H_INCLUDED
#define COMPUTER_H_INCLUDED

//lower fitness indicates that pot is easier.

struct potInfo {
    int ballOrder;
    std::pair<double,double> aimVector;
    double potAngle;
    double distToObject;
    double distToPocket;
    double pocketAngle;
    double fitness;
};

struct _shotInfo {
    double angle;
    double speed;
    double offset;
    double theta;
    int ballOrder;
};

class Computer : public Server
{
    public:
        std::vector<potInfo> potsOn;

        Computer() {};
        double getAngle(double x1, double y1, double x2, double y2);
        bool isBetweenTwoPoints(double x1, double y1, double x2, double y2, double x3, double y3);
        bool isBallIntersectingLine(double x1, double y1, double x2, double y2, double x3, double y3);
        bool isPotBlocked(double x0, double y0, double x1, double y1, int objectBallIndex, double x2, double y2, double x3, double y3);
        void updatePotFitness();
        void getPotOn(int ballOrder);
        void getPotsOn(int ballOrder, bool anyColour);
        _shotInfo getShot(int ballOrder);
};

double Computer::getAngle(double x1, double y1, double x2, double y2)
{
    double d=std::sqrt(x1*x1+y1*y1);
    x1/=d; y1/=d;

    d=std::sqrt(x2*x2+y2*y2);
    x2/=d; y2/=d;

    return std::acos(x1*x2+y1*y2);
}

bool Computer::isBetweenTwoPoints(double x1, double y1, double x2, double y2, double x3, double y3)
{
    //check if 3 is between 1 and 2.

    double dx1=x2-x1; double dy1=y2-y1;
    double dx2=x3-x1; double dy2=y3-y1;
    double dx3=x3-x2; double dy3=y3-y2;

    double a=getAngle(dx1,dy1,dx2,dy2);
    double b=getAngle(-dx1,-dy1,dx3,dy3);

    bool isBetween=false;
    if (a<=pi/2. && b<=pi/2.) {isBetween=true;}

    return isBetween;
}

bool Computer::isBallIntersectingLine(double x1, double y1, double x2, double y2, double x3, double y3)
{
    //1 and 2 define the line.
    //3 is position of ball.

    bool isIntersecting=false;

    bool isBetween=isBetweenTwoPoints(x1,y1,x2,y2,x3,y3);

    if (isBetween)
    {
        isIntersecting=is_ball_blocking(x1,y1,x2,y2,x3,y3);
    }
    else
    {
        //check the end points.
        double d1=std::sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1));
        double d2=std::sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));

        if (d1<=2.*ball_radius || d2<=2.*ball_radius) {isIntersecting=true;}
    }

    return isIntersecting;
}

bool Computer::isPotBlocked(double x0, double y0, double x1, double y1, int objectBallIndex, double x2, double y2, double x3, double y3)
{
    //0 and 1 are start/final positions of cue ball respectively.
    //2 and 3 are start/final positions of object ball respectively.

    bool isBlocked=false;

    for (int i=1;i<22;i++)
    {
        if (serverballs[i]._potted || i==objectBallIndex) {continue;}

        //check if blocking path of cue ball.
        isBlocked=isBallIntersectingLine(x0,y0,x1,y1,serverballs[i]._x,serverballs[i]._y);
        if (isBlocked) {break;}

        //check if blocking path of object ball.
        isBlocked=isBallIntersectingLine(x2,y2,x3,y3,serverballs[i]._x,serverballs[i]._y);
        if (isBlocked) {break;}
    }

    return isBlocked;
}

void Computer::updatePotFitness()
{
    for (int i=0;i<potsOn.size();i++)
    {
        potsOn[i].fitness=(potsOn[i].distToObject+potsOn[i].distToPocket)/100.+potsOn[i].potAngle/(0.5*pi);
    }
}

void Computer::getPotOn(int ballOrder)
{
    double x0=serverballs[0]._x; double y0=serverballs[0]._y;
    double x1, y1; double x2, y2; double x3, y3;

    double dx, dy; double d;

    bool isBlocked;

    x2=serverballs[ballOrder-1]._x; y2=serverballs[ballOrder-1]._y;

    //corner pockets.
    for (int j=0;j<4;j++)
    {
        //set up trajectory.
        bool isBetween=isBetweenTwoPoints(cPottingPoints[j][0],cPottingPoints[j][1],cPottingPoints2[j][0],cPottingPoints2[j][1],x2,y2);
        if (isBetween)
        {
            x3=cPottingPoints2[j][0]; y3=cPottingPoints2[j][1];
        }
        else
        {
            x3=cPottingPoints[j][0]; y3=cPottingPoints[j][1];
        }

        dx=x3-x2; dy=y3-y2;
        d=std::sqrt(dx*dx+dy*dy);
        dx/=d; dy/=d;

        x1=x2-dx*2.*ball_radius; y1=y2-dy*2.*ball_radius;

        //check if overcut.
        d=getAngle(x1-x0,y1-y0,x3-x2,y3-y2);
        if (d>=pi/2.) {continue;}

        //check if blocked.
        isBlocked=isPotBlocked(x0,y0,x1,y1,ballOrder-1,x2,y2,x3,y3);
        if (isBlocked) {continue;}

        //pot is on now.
        potsOn.push_back(potInfo());
        d=std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
        potsOn.back().aimVector=std::pair<double,double>((x1-x0)/d,(y1-y0)/d);
        potsOn.back().ballOrder=ballOrder;
        potsOn.back().distToObject=d;
        d=std::sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
        potsOn.back().distToPocket=d;
        d=getAngle(x1-x0,y1-y0,x3-x2,y3-y2);
        potsOn.back().potAngle=d;

        if (isBetween) {potsOn.back().pocketAngle=0.;}
        else
        {
            potsOn.back().pocketAngle=getAngle(x3-cPottingPoints2[j][0],y3-cPottingPoints2[j][1],x2-x3,y2-y3);
        }
    }

    //middle pockets.
    for (int j=0;j<2;j++)
    {
        //set up trajectory.
        x3=mPottingPoints[j][0]; y3=mPottingPoints[j][1];

        dx=x3-x2; dy=y3-y2;
        d=std::sqrt(dx*dx+dy*dy);
        dx/=d; dy/=d;

        x1=x2-dx*2.*ball_radius; y1=y2-dy*2.*ball_radius;

        //check if angle to pocket too extreme.
        d=getAngle(blue_x-x3,blue_y-y3,x2-x3,y2-y3);
        if (d>=mPocketMaxAngle) {continue;}

        //check if overcut.
        d=getAngle(x1-x0,y1-y0,x3-x2,y3-y2);
        if (d>=pi/2.) {continue;}

        //check if blocked.
        isBlocked=isPotBlocked(x0,y0,x1,y1,ballOrder-1,x2,y2,x3,y3);
        if (isBlocked) {continue;}

        //pot is on now.
        potsOn.push_back(potInfo());
        d=std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
        potsOn.back().aimVector=std::pair<double,double>((x1-x0)/d,(y1-y0)/d);
        potsOn.back().ballOrder=ballOrder;
        potsOn.back().distToObject=d;
        d=std::sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
        potsOn.back().distToPocket=d;
        d=getAngle(x1-x0,y1-y0,x3-x2,y3-y2);
        potsOn.back().potAngle=d;
        potsOn.back().pocketAngle=getAngle(blue_x-x3,blue_y-y3,x2-x3,y2-y3);
    }
}

void Computer::getPotsOn(int ballOrder=8, bool anyColour=true)
{
    potsOn.clear();

    if (ballOrder>=8)
    {
        for (int i=7;i<22;i++)
        {
            getPotOn(i+1);
        }
    }
    else
    {
        if (anyColour)
        {
            for (int i=1;i<7;i++)
            {
                getPotOn(i+1);
            }
        }
        else {getPotOn(ballOrder);}
    }

    updatePotFitness();

    return;
}

_shotInfo Computer::getShot(int ballOrder)
{
    getPotsOn(ballOrder);

    _shotInfo shot;

    if (potsOn.size()==0)
    {
        //play safety shot.
        shot.offset=0.;
        shot.theta=0.;
        shot.speed=45.;
        shot.angle=-0.5*pi;
    }
    else
    {
        //find the easiest pot and play it.
        int index=0;
        double topfit=9999999.;
        for (int i=0;i<potsOn.size();i++)
        {
            if (potsOn[i].fitness<topfit)
            {
                topfit=potsOn[i].fitness;
                index=i;
            }
        }

        //get angle of shot.
        shot.offset=0.;
        shot.theta=0.;
        shot.speed=45.;
        shot.angle=0.5*pi-std::atan2(potsOn[index].aimVector.second,potsOn[index].aimVector.first);
        shot.ballOrder=potsOn[index].ballOrder;
    }

    return shot;
}

#endif // COMPUTER_H_INCLUDED
