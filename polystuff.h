#ifndef POLYSTUFF_H_INCLUDED
#define POLYSTUFF_H_INCLUDED

const double DOUBLE_EPSILON=std::numeric_limits<double>::epsilon();

const int qlimit=100;

template <typename T> int sgn(T val)
{
    return (T(0)<val)-(val<T(0));
}

double nroot(double A,int n)
{
    if (fabs(A)<DOUBLE_EPSILON)
    {
        return double(0.);
    }

    double x_=pow(A,1./n);
    double xi=(1./n)*((n-1)*x_+A/pow(x_,n-1.));
    int oscillate=0;
    int runs=0;
    bool bisection=false;
    bool converged=false;
    double s;
    double u;
    double z;
    while (!converged && !bisection)
    {
        runs+=1;
        if ((pow(x_,n)-A)*(pow(xi,n)-A)<0)
        {
            if (pow(xi,n)-A<0)
            {
                oscillate+=1;
                s=xi;
            }
            else
            {
                u=xi;
            }
        }
        z=(pow(xi,n)-A)/(n*pow(xi,n-1));
        x_=xi;
        //xi=(1./n)*((n-1)*x_+A/pow(x_,n-1.));
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
        double t=u-s;
        while (fabs(t)>fabs(xi)*DOUBLE_EPSILON)
        {
            runs+=1;
            if (pow(xi,n)-A<0)
            {
                s=xi;
            }
            else
            {
                u=xi;
            }
            t=0.5*(u-s);
            xi=s+t;
            if (runs>qlimit)
            {
                break;
            }
        }
    }

    return xi;
}

std::array<double,2> qsolve_quadratic(double a,double b,double c)
{
    std::array<double,2> roots;

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
            roots[0]=sqrt(-1.);
            roots[1]=sqrt(-1.);
        }
        else
        {
            roots[0]=nroot(c,2);
            roots[1]=-nroot(c,2);
        }
        return roots;
    }

    if (pow(b,2)-4*a*c<0)
    {
        //complex roots.
        roots[0]=sqrt(-1.);
        roots[1]=sqrt(-1.);
    }
    else
    {
        double q=-0.5*(b+sgn(b)*sqrt(pow(b,2)-4.*a*c));
        roots[0]=q;
        roots[1]=c/q;
    }
    return roots;
}

std::array<double,3> qsolve_cubic(double a,double b,double c,double d)
{
    std::array<double,3> roots;
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

    double temp=nroot(fabs(a1),2);
    if (temp>k)
    {
        k=temp;
        chosen=1;
    }
    temp=nroot(fabs(a0),3);
    if (temp>k)
    {
        k=temp;
        chosen=0;
    }

    a0=a0/pow(k,3);
    a1=a1/pow(k,2);
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
        roots[0]=sqrt(-1.);
        roots[1]=sqrt(-1.);
        roots[2]=sqrt(-1.);
        return roots;
    }

    //get NR starting point.
    double p;
    double q;
    double r;
    double s;
    double x;
    double xshift=0;

    if (kind==1)
    {
        p=1.09574;
        q=-0.3239;
        r=-0.3239;
        s=0.0957439;
        x=p+q*a1+r*a2+s*a1*a2;
    }
    else if (kind==2)
    {
        p=-1.09574;
        q=0.3239;
        r=-0.3239;
        s=0.0957439;
        x=p+q*a1+r*a2+s*a1*a2;
    }
    else if (kind==3)
    {
        if (a0<-a2/3.+pow(a2,3)*2./27.)
        {
            p=1.14413;
            q=-0.275509;
            r=-0.445578;
            s=-0.0259342;
        }
        else
        {
            p=-1.14413;
            q=-0.275509;
            r=-0.445578;
            s=0.0259342;
        }
        x=p+q*a0+r*a2+s*a0*a2;
    }
    else if (kind==4)
    {
        if (a0>0)
        {
            p=-0.771845;
            q=-0.228155;
        }
        else
        {
            p=-0.771845;
            q=0.228155;
        }
        x=p*a0+q*a0*a2;
    }
    else if (kind==5)
    {
        if (fabs(a1-1./3.)<DOUBLE_EPSILON && fabs(a0+1./27.)<DOUBLE_EPSILON)
        {
            roots[0]=(1./3.)*k;
            roots[1]=(1./3.)*k;
            roots[2]=(1./3.)*k;
            return roots;
        }

        if (-1<=a1 && a1<=1./3.)
        {
            if (a0<-a1/3.+2./27.)
            {
                p=0.878558;
                q=-0.571888;
                r=-0.711154;
                s=-0.322313;
            }
            else
            {
                p=-0.192823;
                q=-0.566324;
                r=0.505734;
                s=-0.264881;
            }
        }
        else if (1./3.<a1 && a1<=1)
        {
            if (a0<-a1/3.+2./27.)
            {
                p=1.19748;
                q=-0.283772;
                r=-0.837476;
                s=-0.356228;
            }
            else
            {
                p=-0.345219;
                q=-0.401231;
                r=0.207216;
                s=-0.00445532;
            }

        }
        x=p+q*a0+r*a1+s*a0*a1;

        if (fabs(a1-1./3.)<=0.01 && fabs(a0+1./27.)<=0.01)
        {
            xshift=-1./3.;
            x=x+xshift;
            a2=0.;
            a1=a1-1./3.;
            a0=a1/3.+a0+1./27.;
            if (fabs(a0)<DOUBLE_EPSILON) {a0=0.;}
        }
    }
    else if (kind==6)
    {
        if (fabs(a1-1./3.)<DOUBLE_EPSILON && fabs(a0-1./27.)<DOUBLE_EPSILON)
        {
            roots[0]=-(1./3.)*k;
            roots[1]=-(1./3.)*k;
            roots[2]=-(1./3.)*k;
            return roots;
        }

        if (-1<=a1 && a1<=1./3.)
        {
            if (a0>a1/3.-2./27.)
            {
                p=-0.878558;
                q=-0.571888;
                r=0.711154;
                s=-0.322313;
            }
            else
            {
                p=0.192823;
                q=-0.566324;
                r=-0.505734;
                s=-0.264881;
            }
        }
        else if (1./3.<a1 && a1<=1)
        {
            if (a0>a1/3.-2./27.)
            {
                p=-1.19748;
                q=-0.283772;
                r=0.837476;
                s=-0.356228;
            }
            else
            {
                p=0.345219;
                q=-0.401231;
                r=-0.207216;
                s=-0.00445532;
            }
        }
        x=p+q*a0+r*a1+s*a0*a1;

        if (fabs(a1-1./3.)<=0.01 && fabs(a0-1./27.)<=0.01)
        {
            xshift=1./3.;
            x=x+xshift;
            a2=0.;
            a1=a1-1./3.;
            a0=a1/3.+a0-1./27.;
            if (fabs(a0)<DOUBLE_EPSILON) {a0=0.;}
        }
    }

    //do a NR to get a real root.
    double x_=x;
    double xi=x_-(pow(x_,3)+a2*pow(x_,2)+a1*x_+a0)/(3*pow(x_,2)+2*a2*x_+a1);
    int oscillate=0;
    int runs=0;
    bool bisection=false;
    bool converged=false;
    double s_;
    double u;
    double z;
    double fx_;
    double fxi;
    while (!converged && !bisection)
    {
        runs+=1;
        fx_=pow(x_,3)+a2*pow(x_,2)+a1*x_+a0;
        fxi=pow(xi,3)+a2*pow(xi,2)+a1*xi+a0;
        if (fx_*fxi<0)
        {
            if (fxi<0)
            {
                oscillate+=1;
                s_=xi;
            }
            else
            {
                u=xi;
            }
        }
        z=(pow(xi,3)+a2*pow(xi,2)+a1*xi+a0)/(3*pow(xi,2)+2*a2*xi+a1);
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
            fxi=pow(xi,3)+a2*pow(xi,2)+a1*xi+a0;
            if (fxi<0)
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

    if (fabs(a1*x)>big)
    {
        i=1;
        big=fabs(a1*x);
    }
    if (fabs(a2*x*x)>big)
    {
        i=2;
        big=fabs(a2*x*x);
    }
    if (fabs(x*x*x)>big)
    {
        i=3;
    }

    x=x*k;
    //get the deflated quadratic from unscaled cubic.
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
    std::array<double,4> roots;

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
        if (qroots[0]<0)
        {
            roots[0]=sqrt(-1.);
            roots[1]=sqrt(-1.);
        }
        else
        {
            roots[0]=nroot(qroots[0],2);
            roots[1]=-nroot(qroots[0],2);
        }
        if (qroots[1]<0)
        {
            roots[2]=sqrt(-1.);
            roots[3]=sqrt(-1.);
        }
        else
        {
            roots[2]=nroot(qroots[1],2);
            roots[3]=-nroot(qroots[1],2);
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

    double temp=nroot(fabs(a2),2);
    if (temp>k)
    {
        k=temp;
        chosen=2;
    }
    temp=nroot(fabs(a1),3);
    if (temp>k)
    {
        k=temp;
        chosen=1;
    }
    temp=nroot(fabs(a0),4);
    if (temp>k)
    {
        k=temp;
        chosen=0;
    }

    a0=a0/pow(k,4);
    a1=a1/pow(k,3);
    a2=a2/pow(k,2);
    a3=a3/k;

    if (chosen==0) {a0=1.*sgn(a0);}
    else if (chosen==1) {a1=1.*sgn(a1);}
    else if (chosen==2) {a2=1.*sgn(a2);}
    else if (chosen==3) {a3=1.*sgn(a3);}

    std::array<double,3> sp=qsolve_cubic(1.,0.75*a3,0.5*a2,0.25*a1);

    double s=sp[0];
    double u=pow(10,308);
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
    if (uroot) {qu=pow(u,4)+a3*pow(u,3)+a2*pow(u,2)+a1*u+a0;}
    else {qu=1.;}
    double qs=pow(s,4)+a3*pow(s,3)+a2*pow(s,2)+a1*s+a0;

    if (qu>0 && qs>0)
    {
        //complex roots only.
        roots[0]=sqrt(-1.);
        roots[1]=sqrt(-1.);
        roots[2]=sqrt(-1.);
        roots[3]=sqrt(-1.);
        return roots;
    }
    else if (qs<0 && qu<0)
    {
        if (qs<qu)
        {
            if (s<0 && a0>0) {x=0;}
            else {x=2.;}
        }
        else
        {
            if (u>0 && a0>0) {x=0;}
            else {x=-2.;}
        }
    }
    else if (qs<0 && qu>=0)
    {
        if (s<-a3/4.)
        {
            if (s>0 && a0>0) {x=0;}
            else {x=-2.;}
        }
        else
        {
            if (s<0 && a0>0) {x=0;}
            else {x=2.;}
        }
    }
    else if (qs>=0 && qu<0)
    {
        if (u<-a3/4.)
        {
            if (u>0 && a0>0) {x=0;}
            else {x=-2.;}
        }
        else
        {
            if (u<0 && a0>0) {x=0;}
            else {x=2.;}
        }
    }

    //do a NR to get a real root.
    double x_=x;
    double xi=x_-(pow(x_,4)+a3*pow(x_,3)+a2*pow(x_,2)+a1*x_+a0)/(4*pow(x_,3)+3*a3*pow(x_,2)+2*a2*x_+a1);
    int oscillate=0;
    int runs=0;
    bool bisection=false;
    bool converged=false;
    double s_;
    double u_;
    double z;
    double fx_;
    double fxi;
    while (!converged && !bisection)
    {
        runs+=1;
        fx_=pow(x_,4)+a3*pow(x_,3)+a2*pow(x_,2)+a1*x_+a0;
        fxi=pow(xi,4)+a3*pow(xi,3)+a2*pow(xi,2)+a1*xi+a0;
        if (fx_*fxi<0)
        {
            if (fxi<0)
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
        z=(pow(x_,4)+a3*pow(x_,3)+a2*pow(x_,2)+a1*x_+a0)/(4*pow(x_,3)+3*a3*pow(x_,2)+2*a2*x_+a1);
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
        double t=u_-s_;
        while (fabs(t)>fabs(xi)*DOUBLE_EPSILON)
        {
            runs+=1;
            fxi=pow(xi,4)+a3*pow(xi,3)+a2*pow(xi,2)+a1*xi+a0;
            if (fxi<0)
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

    if (fabs(a1*x)>big)
    {
        i=1;
        big=fabs(a1*x);
    }
    if (fabs(a2*x*x)>big)
    {
        i=2;
        big=fabs(a2*pow(x,2));
    }
    if (fabs(a3*pow(x,3))>big)
    {
        i=3;
        big=fabs(a3*pow(x,3));
    }
    if (fabs(pow(x,4))>big)
    {
        i=4;
    }

    x=x*k;
    //get the deflated cubic from unscaled quartic.
    double t;
    if (i==4 || i==3)
    {
        u=-e/x;
        t=(u-d)/x;
        s=(t-c)/x;
    }
    else if (i==2)
    {
        u=-e/x;
        t=(u-d)/x;
        s=b+x;
    }
    else if (i==1)
    {
        s=b+x;
        t=c+s*x;
        u=-e/x;
    }
    else if (i==0)
    {
        s=b+x;
        t=c+s*x;
        u=d+t*x;
    }

    std::array<double,3> qroots=qsolve_cubic(1.,s,t,u);

    roots[0]=x;
    roots[1]=qroots[0];
    roots[2]=qroots[1];
    roots[3]=qroots[2];

    return roots;
}

#endif // POLYSTUFF_H_INCLUDED
