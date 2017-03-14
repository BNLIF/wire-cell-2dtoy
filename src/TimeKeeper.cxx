#include "WireCell2dToy/TimeKeeper.h"

#include <sstream>
#include <iostream>		// debugging

using namespace std;
using namespace WireCell;

TimeKeeper::TimeKeeper(const std::string& msg, ptime starting_time)
{
    //cerr << "TimeKeeper starting with \"" << msg << "\"" << endl;
    m_events.push_back(event(starting_time, msg));
}

TimeKeeper::~TimeKeeper()
{
}

std::string TimeKeeper::operator()(std::string msg, ptime now)
{
    //cerr << "TimeKeeper: accepting message: " << msg << endl;
    m_events.push_back(event(now, msg));
    return emit(-1);
}

boost::posix_time::ptime TimeKeeper::start_time() const
{
    return (*this)[0].first;
}
boost::posix_time::ptime TimeKeeper::last_time() const
{
    return (*this)[-1].first;
}
TimeKeeper::deltat TimeKeeper::last_duration() const
{
    return (*this)[-1].first - (*this)[-2].first;
}

boost::posix_time::time_duration TimeKeeper::since(ptime now) const
{
    return now - start_time();
}

TimeKeeper::event TimeKeeper::operator[](int ind) const
{
    while (ind < 0) { ind += m_events.size();}

    return m_events[ind];
}

std::string TimeKeeper::summary() const
{
    stringstream ss;
    for (int ind=0; ind<m_events.size(); ++ind) {
	ss << this->emit(ind) << "\n";
    }
    return ss.str();
}

std::string TimeKeeper::emit(int ind) const
{
    while (ind < 0) { ind += m_events.size();}
    //cerr << "emit " << ind << endl;
    int prev_ind = ind-1;
    if (prev_ind<0) prev_ind=0;
    const event& prev = (*this)[prev_ind];
    const event& evt = (*this)[ind];

    deltat from_start = since(evt.first);
    deltat from_last = evt.first - prev.first;

    stringstream ss;
    ss << "TICK: " << from_start.total_milliseconds() << " ms "
       << "(this: " << from_last.total_milliseconds() << " ms) "
       << evt.second;
    return ss.str();
}
