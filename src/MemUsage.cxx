#include "WireCell2dToy/MemUsage.h"

#include <unistd.h>
#include <iostream>		// debugging
#include <sstream>

using namespace std;
using namespace WireCell;

static double memusage_linux_resident() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double rss = resident * page_size_kb;
    return rss;
}

static double memusage_linux_shared() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double shm = share * page_size_kb;
    return shm;
}

static double memusage_linux_size() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double siz = tSize * page_size_kb;
    return siz;
}


double WireCell::memusage_resident() {
#ifdef __linux__
    return memusage_linux_resident();
#endif
    return -1;
}

double WireCell::memusage_shared()
{
#ifdef __linux__
    return memusage_linux_shared();
#endif
    return -1;
}
double WireCell::memusage_size()
{
#ifdef __linux__
    return memusage_linux_size();
#endif
    return -1;
}


MemUsage::MemUsage(const std::string& msg)
{
    push(msg);
}

MemUsage::~MemUsage()
{
}

MemUsage::memusage MemUsage::current() const
{
    return memusage(memusage_size(), memusage_resident());
}
void MemUsage::push(const std::string& msg, MemUsage::memusage mu)
{
    if (mu.first < 0 && mu.second < 0) {
	mu = current();
    }
    m_events.push_back(event(mu,msg));
}

std::string MemUsage::operator()(std::string msg, MemUsage::memusage mu)
{
    push(msg, mu);
    return emit(-1);
}

	/// Return event by index.
MemUsage::event MemUsage::operator[](int ind) const
{
    while (ind < 0) { ind += m_events.size();}

    return m_events[ind];
}

std::string MemUsage::summary() const
{
    stringstream ss;
    for (int ind=0; ind<m_events.size(); ++ind) {
	ss << this->emit(ind) << "\n";
    }
    return ss.str();
}

std::string MemUsage::emit(int ind) const
{
    while (ind < 0) { ind += m_events.size();}
    int prev_ind = ind-1;
    if (prev_ind<0) prev_ind=0;

    const memusage& prev_mem = (*this)[prev_ind].first;
    const memusage& evt_mem = (*this)[ind].first;
    const string& evt_msg = (*this)[ind].second;

    memusage from_prev(evt_mem.first - prev_mem.first, evt_mem.second - prev_mem.second);

    stringstream ss;
    ss << "MEM: total: size=" << evt_mem.first << "K, res=" << evt_mem.second << "K "
       << "increment: size=" << from_prev.first << "K, res=" << from_prev.second << "K "
       << evt_msg;
    return ss.str();
}
