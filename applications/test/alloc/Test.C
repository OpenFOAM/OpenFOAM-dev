#include <stdlib.h>

class Int
{
    int I;

public:

    Int(){}

    operator int()
    {
        return I;
    }
};


template<class T>
class List : public T
{
    T* v;
    int sz;

public:

    List()
    {
        v = new T[sz=10];
    }

    List(int s)
    {
        v = new T[sz=s];
    }

    ~List()
    {
        delete[] v;
    }

    inline int size() const;

};


template<class T>
inline int List<T>::size() const
{
    return sz;
}


#include <stream.h>

main()
{
    typedef List<Int> intList;

    intList list(10);

    cout << list.size() << "\n";

    return 0;
}
