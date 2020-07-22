/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::CorrectionLimitingMethods::relative

Description
    Correction limiting method based on the relative particle velocity.

    This method limits the velocity correction to that of a rebound with a
    coefficient of restitution \f$e\f$. The relative velocity of the particle
    with respect to the mean value is used to calculate the direction and
    magnitude of the limited velocity.

SourceFiles
    relative.C

\*---------------------------------------------------------------------------*/

#ifndef relative_H
#define relative_H

#include "CorrectionLimitingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace CorrectionLimitingMethods
{

/*---------------------------------------------------------------------------*\
                          Class relative Declaration
\*---------------------------------------------------------------------------*/

class relative
:
    public CorrectionLimitingMethod
{
protected:

    // Protected data

        //- Coefficient of restitution
        scalar e_;


public:

    //- Runtime type information
    TypeName("relative");


    // Constructors

        //- Construct from components
        relative(const dictionary& dict);

        //- Copy constructor
        relative(const relative& cl);

        //- Construct and return a clone
        virtual autoPtr<CorrectionLimitingMethod> clone() const
        {
            return autoPtr<CorrectionLimitingMethod>
            (
                new relative(*this)
            );
        }


    //- Destructor
    virtual ~relative();


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
