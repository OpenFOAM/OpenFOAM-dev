/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

Class
    Foam::functionObjects::cylindrical

Description
    Transforms the specified velocity field into a
    cylindrical polar coordinate system or back to Cartesian.

    Example of function object specification to convert the velocity field U
    into cylindrical polar coordinates before averaging and returning the
    average to Cartesian coordinates:
    \verbatim
    cartesianToCylindrical
    {
        type        cylindrical;
        libs        ("libfieldFunctionObjects.so");

        origin      (0 0 0);
        axis        (0 0 1);

        field       U;

        writeControl    outputTime;
        writeInterval   1;
    }

    #includeFunc fieldAverage(cylindrical(U))

    cylindricalToCartesian
    {
        type        cylindrical;
        libs        ("libfieldFunctionObjects.so");

        origin      (0 0 0);
        axis        (0 0 1);

        field       cylindrical(U)Mean;
        toCartesian true;
        result      UMean;

        writeControl    outputTime;
        writeInterval   1;
    }
    \endverbatim
    This is particularly useful for cases with rotating regions, e.g. mixer
    vessels with NCC.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    cylindrical.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalFunctionObject_H
#define cylindricalFunctionObject_H

#include "fieldExpression.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class cylindrical Declaration
\*---------------------------------------------------------------------------*/

class cylindrical
:
    public fieldExpression
{
    // Private Data

        vector origin_;

        vector axis_;


    // Private Member Functions

        tensor R(const vector& p) const;

        void transform(vectorField& vf, const vectorField& points) const;

        bool toCartesian_;

        //- Calculate the cylindrical field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("cylindrical");


    // Constructors

        //- Construct from Time and dictionary
        cylindrical
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cylindrical();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
