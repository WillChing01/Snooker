#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

//grid squares.
std::array<std::array<std::array<int,7>,39>,73> grid={};
std::array<std::array<int,39>,73> grid_index={};

//matrix.
Eigen::Matrix<double,46,46> M_=Eigen::MatrixXd::Zero(46,46);

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

#endif // COMMON_H_INCLUDED
