#include "scalar.H"
#include "vector.H"
#include "curveTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar distance(const vector& p1, const vector& p2)
{
    return mag(p2 - p1);
}


bool stepForwardsToNextPoint
(
    const vector& o,
    vector& n,
    label& i,
    label& ip1,
    scalar l,
    const curve& Curve
)
{
    label ip1n = ip1-1;
    while (++ip1n < Curve.size() && distance(o, Curve[ip1n]) < l);
    label in = ip1n - 1;

    bool eoc = true;

    if (ip1n < Curve.size() && in >= 0)
    {
        eoc = interpolate(Curve[in], Curve[ip1n], o, n, l);

        i = in;
        ip1 = ip1n;
    }

    return eoc;
}


bool stepBackwardsToNextPoint
(
    const vector& o,
    vector& n,
    label& i,
    label& ip1,
    scalar l,
    const curve& Curve
)
{
    label ip1n = ip1+1;
    while (--ip1n >= 0 && distance(o, Curve[ip1n]) < l);
    label in = ip1n + 1;

    bool eoc = true;

    if (ip1n >= 0 && in < Curve.size())
    {
        eoc = interpolate(Curve[in], Curve[ip1n], o, n, l);

        i = in;
        ip1 = ip1n;
    }

    return eoc;
}


bool interpolate
(
    const vector& p1,
    const vector& p2,
    const vector& o,
    vector& n,
    scalar l
)
{
    vector D = p1 - p2;
    scalar a = magSqr(D);
    scalar b = 2.0*(D&(p2 - o));
    scalar c = magSqr(p2) + (o&(o - 2.0*p2)) - l*l;

    scalar b2m4ac = b*b - 4.0*a*c;

    if (b2m4ac >= 0.0)
    {
        scalar srb2m4ac = sqrt(b2m4ac);

        scalar lambda = (-b - srb2m4ac)/(2.0*a);

        if (lambda > 1.0+curveSmall || lambda < -curveSmall)
        {
            lambda = (-b + srb2m4ac)/(2.0*a);
        }

        if (lambda < 1.0+curveSmall && lambda > -curveSmall)
        {
            n = p2 + lambda*(p1 - p2);

            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        return true;
    }
}



bool XstepForwardsToNextPoint
(
    const vector& o,
    vector& n,
    label& i,
    label& ip1,
    scalar l,
    const curve& Curve
)
{
    label ip1n = ip1-1;
    while (++ip1n < Curve.size() && mag(o.x() - Curve[ip1n].x()) < l);
    label in = ip1n - 1;

    bool eoc = true;

    if (ip1n < Curve.size() && in >= 0)
    {
        eoc = Xinterpolate(Curve[in], Curve[ip1n], o, n, l);

        i = in;
        ip1 = ip1n;
    }

    return eoc;
}



bool Xinterpolate
(
    const vector& p1,
    const vector& p2,
    const vector& o,
    vector& n,
    scalar l
)
{
    n.x() = o.x() + l;

    if (p2.x() < o.x() + l - curveSmall && p2.x() > o.x() - l + curveSmall)
    {
        return true;
    }

    if (p2.x() < o.x() + l)
    {
        n.x() = o.x() - l;
    }

    vector D = p2 - p1;
    scalar lambda = (n.x() - p1.x())/D.x();
    n.y() = p1.y() + lambda*D.y();
    return false;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
