#ifndef WIRECELLUTIL_TIMEKEEPER
#define WIRECELLUTIL_TIMEKEEPER

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <vector>

namespace WireCell {

    /** A helper class to give some time keeping.
     *
     * Use like
     *
     *   TimeKeeper tk("starting");
     *   ...
     *   cerr << tk("starting long calculation....") << endl;
     *   do_long_calculation();
     *   cerr << tk("...done") << endl;;
     *   ...
     *   cerr << tk.summary() << endl;
     */
    class TimeKeeper {
    public:
	typedef boost::posix_time::ptime ptime;
	typedef boost::posix_time::time_duration deltat;
	typedef std::pair<ptime, std::string> event;

	TimeKeeper(const std::string& msg = "start",
		   ptime starting_time = boost::posix_time::microsec_clock::local_time());
	~TimeKeeper();

	/// Return the time at which this time keeper was started.
	ptime start_time() const;

	/// Return the time of the last event.
	ptime last_time() const;

	/// Return the duration between the last two events.
	deltat last_duration() const;

	/// Return the time duration between "now" and the start time.
	deltat since(ptime now = boost::posix_time::microsec_clock::local_time()) const;

	/// Record an event.
	std::string operator()(
	    std::string msg = "<tick>",
	    ptime now = boost::posix_time::microsec_clock::local_time());
	
	/// Return summary up to now.
	std::string summary() const;

	/// Return event by index.
	event operator[](int ind) const;


    private:
	/// Emit a formatted message for the given event index.
	std::string emit(int ind) const;


	std::vector< event > m_events;
    };
}
#endif
