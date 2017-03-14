// http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
// http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c

#ifndef WIRECELL_MEMUSAGE
#define WIRECELL_MEMUSAGE

#include <fstream>
#include <vector>
#include <string>


namespace WireCell {

    double memusage_resident();
    double memusage_shared();
    double memusage_size();

    
    class MemUsage {
    public:	
	typedef std::pair<double, double> memusage;
	typedef std::pair<memusage, std::string> event;
	MemUsage(const std::string& msg = "start");
	~MemUsage();

	void push(const std::string& msg, MemUsage::memusage mu = memusage(-1,-1));

	/// Record an event.
	std::string operator()(
	    std::string msg = "<tick>",
	    MemUsage::memusage mu = memusage(-1,-1));
	
	/// Return summary up to now.
	std::string summary() const;

	/// Return event by index.
	event operator[](int ind) const;
	
	memusage current() const;

    private:
	/// Emit a formatted message for the given event index.
	std::string emit(int ind) const;


	std::vector< event > m_events;

    };

}


#endif
