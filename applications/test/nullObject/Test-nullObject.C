#include "nullObject.H"
#include "IOstreams.H"

using namespace Foam;

class SimpleClass
{

public:

    //- Null constructor
    SimpleClass()
    {}
};


int main()
{
    // Test pointer and reference to a class

    SimpleClass* ptrToClass = new SimpleClass;
    SimpleClass& refToClass(*ptrToClass);

    if (notNull(ptrToClass))
    {
        Info<< "Pass: ptrToClass is not null" << endl;
    }
    else
    {
        Info<< "FAIL: refToClass is null" << endl;
    }

    if (notNull(refToClass))
    {
        Info<< "Pass: refToClass is not null" << endl;
    }
    else
    {
        Info<< "FAIL: refToClass is null" << endl;
    }


    // Test pointer and reference to the nullObject

    const SimpleClass* ptrToNull(NullObjectPtr<SimpleClass>());
    const SimpleClass& refToNull(*ptrToNull);

    if (isNull(ptrToNull))
    {
        Info<< "Pass: ptrToNull is null" << endl;
    }
    else
    {
        Info<< "FAIL: ptrToNull is not null" << endl;
    }

    if (isNull(refToNull))
    {
        Info<< "Pass: refToNull is null" << endl;
    }
    else
    {
        Info<< "FAIL: refToNull is not null" << endl;
    }

    // Clean-up
    delete ptrToClass;

    return 0;
}
