#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

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

#endif // COMMON_H_INCLUDED
