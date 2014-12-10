/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::fv::temperatureLimitsConstraint

Description
    Constraint for temperature to apply limits between minimum and maximum
    values

    Constraint described by:

        temperatureLimitsConstraintCoeffs
        {
            minimum     200;
            maximum     500;
        }


SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureLimitsConstraint_H
#define temperatureLimitsConstraint_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class temperatureLimitsConstraint Declaration
\*---------------------------------------------------------------------------*/

class temperatureLimitsConstraint
:
    public option
{

protected:

    // Protected data

        //- Minimum temperature limit [K]
        scalar Tmin_;

        //- Maximum temperature limit [K]
        scalar Tmax_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        temperatureLimitsConstraint(const temperatureLimitsConstraint&);

        //- Disallow default bitwise assignment
        void operator=(const temperatureLimitsConstraint&);


public:

    //- Runtime type information
    TypeName("temperatureLimitsConstraint");


    // Constructors

        //- Construct from components
        temperatureLimitsConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~temperatureLimitsConstraint()
    {}


    // Member Functions

        virtual bool alwaysApply() const;


        // Evaluate

            //- Correct the energy field
            virtual void correct(volScalarField& he);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
