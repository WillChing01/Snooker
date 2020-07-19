#ifndef POLYSTUFF_H_INCLUDED
#define POLYSTUFF_H_INCLUDED

#include <Eigen/Eigenvalues>
#include <math.h>

const double DOUBLE_EPSILON=std::numeric_limits<double>::epsilon();

const int qlimit=100;

template <typename T> int sgn(T val)
{
    return (T(0)<val)-(val<T(0));
}

std::array<double,8> qsolve_octic(std::array<double,8> a)
{
    //assumes it is a monic octic.
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(8,8);

    A(0,7)=-a[0];
    for (int i=1;i<8;i++)
    {
        A(i,7)=-a[i];
        A(i,i-1)=1.;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(A,false);

    std::array<double,8> roots;
    for (int i=0;i<8;i++)
    {
        if (fabs(es.eigenvalues()[i].imag())<DOUBLE_EPSILON)
        {
            //real root.
            roots[i]=es.eigenvalues()[i].real();
        }
        else
        {
            roots[i]=sqrt(-1.);
        }
    }
    return roots;
}

std::array<double,2> qsolve_quadratic(double a,double b,double c)
{
    std::array<double,2> roots;
    static double nan=sqrt(-1.);

    //check if actually a quadratic.
    if (fabs(a)<DOUBLE_EPSILON)
    {
        //linear?
        if (fabs(b)<DOUBLE_EPSILON)
        {
            roots[0]=nan;
        }
        else
        {
            roots[0]=-c/b;
        }
        roots[1]=nan;
        return roots;
    }

    //rescale.
    b=b/a;
    c=c/a;
    a=1.;

    //special cases.
    if (fabs(b)<DOUBLE_EPSILON && fabs(c)<DOUBLE_EPSILON)
    {
        roots[0]=0.;
        roots[1]=0.;
        return roots;
    }
    if (fabs(c)<DOUBLE_EPSILON)
    {
        roots[0]=0.;
        roots[1]=-b;
        return roots;
    }
    if (fabs(b)<DOUBLE_EPSILON)
    {
        if (c>0)
        {
            roots[0]=nan;
            roots[1]=nan;
        }
        else
        {
            roots[0]=sqrt(-c);
            roots[1]=-sqrt(-c);
        }
        return roots;
    }

    double thing=pow(b,2.)-4.*a*c;

    if (thing<0)
    {
        //complex roots.
        roots[0]=nan;
        roots[1]=nan;
    }
    else
    {
        double q=-0.5*(b+sgn(b)*sqrt(thing));
        roots[0]=q;
        roots[1]=c/q;

        //do a NR to get a real root.
        double x_;
        double xi;
        int oscillate=0;
        int runs=0;
        bool bisection=false;
        bool converged=false;
        double s_=0.;
        double u=0.;
        double z;
        double fx_;
        double fxi;
        for (int i=0;i<2;i++)
        {
            x_=roots[i];
            xi=x_-(pow(x_,2.)+b*x_+c)/(2.*x_+b);
            oscillate=0;
            runs=0;
            bisection=false;
            converged=false;
            while (!converged && !bisection)
            {
                runs+=1;
                fx_=pow(x_,2.)+b*x_+c;
                fxi=pow(xi,2.)+b*xi+c;
                if (fx_*fxi<0.)
                {
                    if (fxi<0.)
                    {
                        oscillate+=1;
                        s_=xi;
                    }
                    else
                    {
                        u=xi;
                    }
                }
                z=(pow(xi,2.)+b*xi+c)/(2.*xi+b);
                x_=xi;
                xi=x_-z;
                if (oscillate>2)
                {
                    bisection=true;
                }
                if (fabs(z)<fabs(xi)*DOUBLE_EPSILON)
                {
                    converged=true;
                }
                if (runs>qlimit)
                {
                    break;
                }
            }

            if (bisection)
            {
                runs=0;
                double t=u-s_;
                while (fabs(t)>fabs(xi)*DOUBLE_EPSILON)
                {
                    runs+=1;
                    fxi=pow(xi,2.)+b*xi+c;
                    if (fxi<0.)
                    {
                        s_=xi;
                    }
                    else
                    {
                        u=xi;
                    }
                    t=0.5*(u-s_);
                    xi=s_+t;
                    if (runs>qlimit)
                    {
                        break;
                    }
                }
            }
            roots[i]=xi;
        }
    }
    return roots;
}

std::array<double,3> qsolve_cubic(double a,double b,double c,double d)
{
    std::array<double,3> roots;
    static double third=1./3.;
    static double two27=2./27.;
    static double one27=1./27.;
    static double nan=sqrt(-1.);

    static double p1=1.09574;
    static double q1=-0.3239;
    static double r1=-0.3239;
    static double s1=0.0957439;

    static double p2=-1.09574;
    static double q2=0.3239;
    static double r2=-0.3239;
    static double s2=0.0957439;

    static double p3=1.14413;
    static double q3=-0.275509;
    static double r3=-0.445578;
    static double s3=-0.0259342;

    static double p4=-0.771845;
    static double q4=-0.228155;

    static double p51=0.878558;
    static double q51=-0.571888;
    static double r51=-0.711154;
    static double s51=-0.322313;

    static double p52=-0.192823;
    static double q52=-0.566324;
    static double r52=0.505734;
    static double s52=-0.264881;

    static double p53=1.19748;
    static double q53=-0.283772;
    static double r53=-0.837476;
    static double s53=-0.356228;

    static double p54=-0.345219;
    static double q54=-0.401231;
    static double r54=0.207216;
    static double s54=-0.00445532;

    //check if actually a cubic.
    if (fabs(a)<DOUBLE_EPSILON)
    {
        //a quadratic?
        std::array<double,2> qroots=qsolve_quadratic(b,c,d);
        roots[0]=qroots[0];
        roots[1]=qroots[1];
        roots[2]=nan;
        return roots;
    }

    //turn into monic.
    b=b/a;
    c=c/a;
    d=d/a;
    a=1.;

    //special cases.
    if (fabs(b)<DOUBLE_EPSILON && fabs(c)<DOUBLE_EPSILON && fabs(d)<DOUBLE_EPSILON)
    {
        roots[0]=0.;
        roots[1]=0.;
        roots[2]=0.;
        return roots;
    }
    if (fabs(c)<DOUBLE_EPSILON && fabs(d)<DOUBLE_EPSILON)
    {
        roots[0]=0.;
        roots[1]=0.;
        roots[2]=-b;
        return roots;
    }
    if (fabs(d)<DOUBLE_EPSILON)
    {
        std::array<double,2> qroots=qsolve_quadratic(1.,b,c);
        roots[0]=0.;
        roots[1]=qroots[0];
        roots[2]=qroots[1];
        return roots;
    }

    double a0=d;
    double a1=c;
    double a2=b;

    //perform rescaling.
    double k=fabs(a2);
    int chosen=2;

    double temp=sqrt(fabs(a1));
    if (temp>k)
    {
        k=temp;
        chosen=1;
    }
    temp=cbrt(fabs(a0));
    if (temp>k)
    {
        k=temp;
        chosen=0;
    }

    a0=a0/pow(k,3.);
    a1=a1/pow(k,2.);
    a2=a2/k;

    if (chosen==0) {a0=1.*sgn(a0);}
    else if (chosen==1) {a1=1.*sgn(a1);}
    else if (chosen==2) {a2=1.*sgn(a2);}

    int kind=0;
    //assign a class. can be more than one class.
    if (fabs(a0+1.)<DOUBLE_EPSILON && kind==0) {kind=1;}
    else if (fabs(a0-1.)<DOUBLE_EPSILON && kind==0) {kind=2;}
    if (fabs(a1+1.)<DOUBLE_EPSILON && kind==0) {kind=3;}
    else if (fabs(a1-1.)<DOUBLE_EPSILON && kind==0) {kind=4;}
    if (fabs(a2+1.)<DOUBLE_EPSILON && kind==0) {kind=5;}
    else if (fabs(a2-1.)<DOUBLE_EPSILON && kind==0) {kind=6;}

    if (kind==0)
    {
        std::cout << "ERROR. COULD NOT ASSIGN A CLASS TO CUBIC." << std::endl;
        roots[0]=nan;
        roots[1]=nan;
        roots[2]=nan;
        return roots;
    }

    //get NR starting point.
    double x=0.;
    double xshift=0;

    if (kind==1)
    {
        x=p1+q1*a1+r1*a2+s1*a1*a2;
    }
    else if (kind==2)
    {
        x=p2+q2*a1+r2*a2+s2*a1*a2;
    }
    else if (kind==3)
    {
        if (a0<-a2*third+pow(a2,3.)*two27)
        {
            x=p3+q3*a0+r3*a2+s3*a0*a2;
        }
        else
        {
            x=-p3+q3*a0+r3*a2-s3*a0*a2;
        }
    }
    else if (kind==4)
    {
        if (a0>0.)
        {
            x=p4*a0+q4*a0*a2;
        }
        else
        {
            x=p4*a0-q4*a0*a2;
        }
    }
    else if (kind==5)
    {
        if (fabs(a1-third)<DOUBLE_EPSILON && fabs(a0+one27)<DOUBLE_EPSILON)
        {
            roots[0]=third*k;
            roots[1]=third*k;
            roots[2]=third*k;
            return roots;
        }

        if (-1.<=a1 && a1<=third)
        {
            if (a0<-a1*third+two27)
            {
                x=p51+q51*a0+r51*a1+s51*a0*a1;
            }
            else
            {
                x=p52+q52*a0+r52*a1+s52*a0*a1;
            }
        }
        else if (third<a1 && a1<=1.)
        {
            if (a0<-a1*third+two27)
            {
                x=p53+q53*a0+r53*a1+s53*a0*a1;
            }
            else
            {
                x=p54+q54*a0+r54*a1+s54*a0*a1;
            }

        }

        if (fabs(a1-third)<=0.01 && fabs(a0+one27)<=0.01)
        {
            xshift=-third;
            x=x+xshift;
            a2=0.;
            a1=a1-third;
            a0=a1*third+a0+one27;
            if (fabs(a0)<DOUBLE_EPSILON) {a0=0.;}
        }
    }
    else if (kind==6)
    {
        if (fabs(a1-third)<DOUBLE_EPSILON && fabs(a0-one27)<DOUBLE_EPSILON)
        {
            roots[0]=-third*k;
            roots[1]=-third*k;
            roots[2]=-third*k;
            return roots;
        }

        if (-1.<=a1 && a1<=third)
        {
            if (a0>a1*third-two27)
            {
                x=-p51+q51*a0-r51*a1+s51*a0*a1;
            }
            else
            {
                x=-p52+q52*a0-r52*a1+s52*a0*a1;
            }
        }
        else if (third<a1 && a1<=1.)
        {
            if (a0>a1*third-two27)
            {
                x=-p53+q53*a0-r53*a1+s53*a0*a1;
            }
            else
            {
                x=-p54+q54*a0-r54*a1+s54*a0*a1;
            }
        }

        if (fabs(a1-third)<=0.01 && fabs(a0-one27)<=0.01)
        {
            xshift=third;
            x=x+xshift;
            a2=0.;
            a1=a1-third;
            a0=a1*third+a0-one27;
            if (fabs(a0)<DOUBLE_EPSILON) {a0=0.;}
        }
    }

    //do a NR to get a real root.
    double x_=x;
    double xi=x_-(pow(x_,3.)+a2*pow(x_,2.)+a1*x_+a0)/(3.*pow(x_,2.)+2.*a2*x_+a1);
    int oscillate=0;
    int runs=0;
    bool bisection=false;
    bool converged=false;
    double s_=0.;
    double u=0.;
    double z;
    double fx_;
    double fxi;
    while (!converged && !bisection)
    {
        runs+=1;
        fx_=pow(x_,3.)+a2*pow(x_,2.)+a1*x_+a0;
        fxi=pow(xi,3.)+a2*pow(xi,2.)+a1*xi+a0;
        if (fx_*fxi<0.)
        {
            if (fxi<0.)
            {
                oscillate+=1;
                s_=xi;
            }
            else
            {
                u=xi;
            }
        }
        z=(pow(xi,3.)+a2*pow(xi,2.)+a1*xi+a0)/(3.*pow(xi,2.)+2.*a2*xi+a1);
        x_=xi;
        xi=x_-z;
        if (oscillate>2)
        {
            bisection=true;
        }
        if (fabs(z)<fabs(xi)*DOUBLE_EPSILON)
        {
            converged=true;
        }
        if (runs>qlimit)
        {
            break;
        }
    }

    if (bisection)
    {
        runs=0;
        double t=u-s_;
        while (fabs(t)>fabs(xi)*DOUBLE_EPSILON)
        {
            runs+=1;
            fxi=pow(xi,3.)+a2*pow(xi,2.)+a1*xi+a0;
            if (fxi<0.)
            {
                s_=xi;
            }
            else
            {
                u=xi;
            }
            t=0.5*(u-s_);
            xi=s_+t;
            if (runs>qlimit)
            {
                break;
            }
        }
    }

    x=xi-xshift;

    int i=0;
    double big=fabs(a0);

    temp=fabs(a1*x);
    if (temp>big)
    {
        i=1;
        big=temp;
    }
    temp=fabs(a2*x*x);
    if (temp>big)
    {
        i=2;
        big=temp;
    }
    if (fabs(x*x*x)>big)
    {
        i=3;
    }

    x=x*k;
    //get the deflated quadratic from unscaled cubic.
    double q=0.;
    double p=0.;
    if (i==3 || i==2)
    {
        q=-d/x;
        p=(q-c)/x;
    }
    else if (i==1)
    {
        p=b+x;
        q=-d/x;
    }
    else if (i==0)
    {
        p=b+x;
        q=c+p*x;
    }

    std::array<double,2> qroots=qsolve_quadratic(1.,p,q);
    roots[0]=x;
    roots[1]=qroots[0];
    roots[2]=qroots[1];

    return roots;
}

std::array<double,4> qsolve_quartic(double a,double b,double c,double d,double e)
{
    std::array<double,4> roots={};
    static double nan=sqrt(-1.);

    bool speak=false;
    if (fabs(a-2.09)<0.01 && fabs(b+66.23)<0.01 && fabs(c-597.90)<0.01 && fabs(d+1186.32)<0.01 && fabs(e-669.24)<0.01)
    {
        speak=true;
    }

    //check if actually a quartic.
    if (fabs(a)<DOUBLE_EPSILON)
    {
        //a cubic?
        std::array<double,3> qroots=qsolve_cubic(b,c,d,e);
        roots[0]=qroots[0];
        roots[1]=qroots[1];
        roots[2]=qroots[2];
        roots[3]=nan;
        return roots;
    }

    //turn into monic.
    b=b/a;
    c=c/a;
    d=d/a;
    e=e/a;
    a=1.;

    //special cases.
    if (fabs(e)<DOUBLE_EPSILON)
    {
        std::array<double,3> qroots=qsolve_cubic(1.,b,c,d);
        roots[0]=0.;
        roots[1]=qroots[0];
        roots[2]=qroots[1];
        roots[3]=qroots[2];
        return roots;
    }
    if (fabs(b)<DOUBLE_EPSILON && fabs(d)<DOUBLE_EPSILON)
    {
        std::array<double,2> qroots=qsolve_quadratic(1.,c,e);
        if (qroots[0]<0.)
        {
            roots[0]=nan;
            roots[1]=nan;
        }
        else
        {
            roots[0]=sqrt(qroots[0]);
            roots[1]=-sqrt(qroots[0]);
        }
        if (qroots[1]<0.)
        {
            roots[2]=nan;
            roots[3]=nan;
        }
        else
        {
            roots[2]=sqrt(qroots[1]);
            roots[3]=-sqrt(qroots[1]);
        }
        return roots;
    }

    double a0=e;
    double a1=d;
    double a2=c;
    double a3=b;

    //perform rescaling.
    double k=fabs(a3);
    int chosen=3;

    double temp=sqrt(fabs(a2));
    if (temp>k)
    {
        k=temp;
        chosen=2;
    }
    temp=cbrt(fabs(a1));
    if (temp>k)
    {
        k=temp;
        chosen=1;
    }
    temp=pow(fabs(a0),0.25);
    if (temp>k)
    {
        k=temp;
        chosen=0;
    }

    a0=a0/pow(k,4.);
    a1=a1/pow(k,3.);
    a2=a2/pow(k,2.);
    a3=a3/k;

    if (chosen==0) {a0=1.*sgn(a0);}
    else if (chosen==1) {a1=1.*sgn(a1);}
    else if (chosen==2) {a2=1.*sgn(a2);}
    else if (chosen==3) {a3=1.*sgn(a3);}

    std::array<double,3> sp=qsolve_cubic(1.,0.75*a3,0.5*a2,0.25*a1);
//
//    if (speak)
//    {
//        std::cout << "TP Coefficients:" << std::endl;
//        std::cout << "b: " << 0.75*a3 << std::endl;
//        std::cout << "c: " << 0.5*a2 << std::endl;
//        std::cout << "d: " << 0.25*a1 << std::endl;
//        std::cout << "Turning points:" << std::endl;
//        std::cout << "Tp: " << sp[0] << std::endl;
//        std::cout << "Tp: " << sp[1] << std::endl;
//        std::cout << "Tp: " << sp[2] << std::endl;
//    }

    double s=sp[0];
    double u=pow(10,308.);
    bool uroot=false;

    for (int i=1;i<3;i++)
    {
        if (sp[i]!=sp[i])
        {
            //nan.
            continue;
        }
        if (sp[i]<u)
        {
            u=sp[i];
            uroot=true;
        }
    }

    if (u>s && uroot)
    {
        double temp=u;
        u=s;
        s=temp;
    }

    double x;
    double qu;
    if (uroot) {qu=pow(u,4.)+a3*pow(u,3.)+a2*pow(u,2.)+a1*u+a0;}
    else {qu=1.;}
    double qs=pow(s,4.)+a3*pow(s,3.)+a2*pow(s,2.)+a1*s+a0;

    if (qu>0. && qs>0.)
    {
        //complex roots only.
        roots[0]=nan;
        roots[1]=nan;
        roots[2]=nan;
        roots[3]=nan;
        return roots;
    }
    else if (qs<0. && qu<0.)
    {
        if (qs<qu)
        {
            if (s<0. && a0>0.) {x=0.;}
            else {x=2.;}
        }
        else
        {
            if (u>0. && a0>0.) {x=0.;}
            else {x=-2.;}
        }
    }
    else if (qs<0. && qu>=0.)
    {
        if (s<-a3/4.)
        {
            if (s>0. && a0>0.) {x=0.;}
            else {x=-2.;}
        }
        else
        {
            if (s<0. && a0>0.) {x=0.;}
            else {x=2.;}
        }
    }
    else if (qs>=0. && qu<0.)
    {
        if (u<-a3/4.)
        {
            if (u>0. && a0>0.) {x=0.;}
            else {x=-2.;}
        }
        else
        {
            if (u<0. && a0>0.) {x=0.;}
            else {x=2.;}
        }
    }

    //do a NR to get a real root.
    double x_=x;
    double xi=x_-(pow(x_,4.)+a3*pow(x_,3.)+a2*pow(x_,2.)+a1*x_+a0)/(4.*pow(x_,3.)+3.*a3*pow(x_,2.)+2.*a2*x_+a1);
    int oscillate=0;
    int runs=0;
    bool bisection=false;
    bool converged=false;
    double s_=0.;
    double u_=0.;
    double z;
    double fx_;
    double fxi;
    while (!converged && !bisection)
    {
//        if (speak)
//        {
//            std::cout << x_ << std::endl;
//        }
        runs+=1;
        fx_=pow(x_,4.)+a3*pow(x_,3.)+a2*pow(x_,2.)+a1*x_+a0;
        fxi=pow(xi,4.)+a3*pow(xi,3.)+a2*pow(xi,2.)+a1*xi+a0;
        if (fx_*fxi<0.)
        {
            if (fxi<0.)
            {
                oscillate+=1;
                s_=xi;
            }
            else
            {
                u_=xi;
            }
        }
        x_=xi;
        z=(pow(x_,4.)+a3*pow(x_,3.)+a2*pow(x_,2.)+a1*x_+a0)/(4.*pow(x_,3.)+3.*a3*pow(x_,2.)+2.*a2*x_+a1);
        xi=x_-z;
        if (oscillate>2)
        {
            bisection=true;
        }
        if (fabs(z)<fabs(xi)*DOUBLE_EPSILON)
        {
            converged=true;
        }
        if (runs>qlimit)
        {
            break;
        }
    }

    if (bisection)
    {
//        if (speak)
//        {
//            std::cout << "Bisection!" << std::endl;
//        }
        runs=0;
        double t=u_-s_;
        while (fabs(t)>fabs(xi)*DOUBLE_EPSILON)
        {
            runs+=1;
            fxi=pow(xi,4.)+a3*pow(xi,3.)+a2*pow(xi,2.)+a1*xi+a0;
            if (fxi<0.)
            {
                s_=xi;
            }
            else
            {
                u_=xi;
            }
            t=0.5*(u_-s_);
            xi=s_+t;
            if (runs>qlimit)
            {
                break;
            }
        }
    }

    x=xi;

    int i=0;
    double big=fabs(a0);

    temp=fabs(a1*x);
    if (temp>big)
    {
        i=1;
        big=temp;
    }
    temp=fabs(a2*x*x);
    if (temp>big)
    {
        i=2;
        big=temp;
    }
    temp=fabs(a3*pow(x,3.));
    if (temp>big)
    {
        i=3;
        big=temp;
    }
    if (fabs(pow(x,4.))>big)
    {
        i=4;
    }

    x=x*k;
    //get the deflated cubic from unscaled quartic.
    double t9;
    double s9;
    double u9;
    if (i==4 || i==3)
    {
        u9=-e/x;
        t9=(u9-d)/x;
        s9=(t9-c)/x;
    }
    else if (i==2)
    {
        u9=-e/x;
        t9=(u9-d)/x;
        s9=b+x;
    }
    else if (i==1)
    {
        s9=b+x;
        t9=c+s9*x;
        u9=-e/x;
    }
    else if (i==0)
    {
        s9=b+x;
        t9=c+s9*x;
        u9=d+t9*x;
    }

//    if (speak)
//    {
//        std::cout << s9 << " : " << t9 << " : " << u9 << std::endl;
//    }

    std::array<double,3> qroots=qsolve_cubic(1.,s9,t9,u9);

    roots[0]=x;
    roots[1]=qroots[0];
    roots[2]=qroots[1];
    roots[3]=qroots[2];

//    if (speak)
//    {
//        std::cout << x << std::endl;
//    }

    return roots;
}

#endif // POLYSTUFF_H_INCLUDED
