/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::CorrectionLimitingMethods::absolute

Description
    Correction limiting method based on the absolute particle velocity.

    This method that limits the velocity correction to that of a rebound with a
    coefficient of restitution \f$e\f$. The absolute velocity of the particle
    is used when calculating the magnitude of the limited correction.
    The direction is calculated using the relative velocity.

SourceFiles
    absolute.C

\*---------------------------------------------------------------------------*/

#ifndef absolute_H
#define absolute_H

#include "CorrectionLimitingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CorrectionLimitingMethods
{

/*---------------------------------------------------------------------------*\
                          Class absolute Declaration
\*---------------------------------------------------------------------------*/

class absolute
:
    public CorrectionLimitingMethod
{
protected:

    // Protected data

        //- Coefficient of restitution
        scalar e_;


public:

    //- Runtime type information
    TypeName("absolute");


    // Constructors

        //- Construct from components
        absolute(const dictionary& dict);

        //- Construct as copy
        absolute(const absolute& cl);

        //- Construct and return a clone
        virtual autoPtr<CorrectionLimitingMethod> clone() const
        {
            return autoPtr<CorrectionLimitingMethod>
            (
                new absolute(*this)
            );
        }


    //- Destructor
    virtual ~absolute();


    // Member Functions

        //- Return the limited velocity
        virtual vector limitedVelocity
        (
            const vector uP,
            const vector dU,
            const vector uMean
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CorrectionLimitingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
