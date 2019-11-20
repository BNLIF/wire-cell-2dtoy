#include "WCPData/Singleton.h"
//#include "WCPUtil/Testing.h"

#include <boost/assert.hpp>

#define Assert BOOST_ASSERT
#define AssertMsg BOOST_ASSERT_MSG 

#include <iostream>
using namespace std;

class Foo {
public:
    Foo() { cout << "Foo() at " << (void*)this << endl; }
    virtual ~Foo() { cout << "~Foo() at " << (void*)this << endl; }
    virtual void chirp() { cout << "Foo::chirp() at " << (void*)this << endl; }
};

typedef WCP::Singleton<Foo> OnlyFoo;

int main()
{
    cout << "First time:" << endl;
    Foo* foo1 = &OnlyFoo::Instance();
    Foo* foo2 = &OnlyFoo::Instance();

    Assert(foo1 == foo2);

    OnlyFoo::Instance().chirp();
    cout << "Second time:" << endl;
    OnlyFoo::Instance().chirp();

    Foo* foo3 = &WCP::Singleton<Foo>::Instance();
    Assert(foo3 == foo1);

}
