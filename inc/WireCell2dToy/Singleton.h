#ifndef WIRECELL_SINGLETON
#define WIRECELL_SINGLETON

namespace WireCell {

    template <class T>
    class Singleton
    {
    public:
	static T& Instance() {
	    static T instance;
	    return instance;
	}

    private:
	Singleton(){}
	~Singleton(){}
	Singleton(Singleton const&){}
	Singleton& operator=(Singleton const&){}
    };
}

#endif
