/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "curve.H"
//#include "curveTools.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// construct as interpolation
/*
curve::curve(const curve& Curve, const label nFacets)
:
    Name("Interpolated" + Curve.Name),
    Style(Curve.Style),
    X(2*nFacets),
    Y(2*nFacets)
{
    // Calculate curve length
    scalar curveLength=0;
    register label i;
    for (i=0; i<Curve.size()-1; i++)
    {
        curveLength += distance(Curve[i], Curve[i+1]);
    }

    scalar stepLength = curveLength/nFacets;
    label nPoints = 0;
    label previous=0, next=1;
    bool endOfCurve;
    vector presentPoint=Curve[0], nextPoint;

    do
    {
        endOfCurve =
        stepForwardsToNextPoint
        (
            presentPoint,
            nextPoint,
            previous,
            next,
            stepLength,
            Curve
        );

        if (!endOfCurve)
        {
            if (nPoints >= size()-1)
            {
                setSize(label(1.5*size()));
            }

            presentPoint = nextPoint;

            x()[nPoints] = nextPoint.x();
            y()[nPoints] = nextPoint.y();

            nPoints++;
        }

    } while (!endOfCurve);

    setSize(nPoints);
}
*/


// construct given name, style and size
curve::curve
(
    const string& name,
    const curveStyle& style,
    const label l
)
:
    scalarField(l, 0.0),
    name_(name),
    style_(style)
{}


// construct from the bits
curve::curve
(
    const string& name,
    const curveStyle& style,
    const scalarField& y
)
:
    scalarField(y),
    name_(name),
    style_(style)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// Gradient operation
/*
curve grad(const curve& Curve)
{
    curve gradCurve(Curve);

    register label i;
    for (i=1; i<Curve.size()-1; i++)
    {
        scalar deltaIm1 = Curve[i].x() - Curve[i-1].x();
        scalar deltaI = Curve[i+1].x() - Curve[i].x();

        scalar deltaAv = 1.0/deltaIm1 + 1.0/deltaI;

        gradCurve.y()[i] =
            (
                (Curve[i+1].y() - Curve[i].y())/sqr(deltaI)
              + (Curve[i].y() - Curve[i-1].y())/sqr(deltaIm1)
            )/deltaAv;
    }

    gradCurve.y()[0] =
        (Curve[1].y() - Curve[0].y())/(Curve[1].x() - Curve[0].x());

    label n = Curve.size()-1;

    gradCurve.y()[n] =
        (Curve[n].y() - Curve[n-1].y())/(Curve[n].x() - Curve[n-1].x());

    return gradCurve;
}
*/


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const curve& c)
{
    os  << nl
        << c.name_ << nl
        << c.style_ << nl
        << static_cast<const scalarField&>(c);

    os.check("Ostream& operator>>(Ostream&, const curve&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
