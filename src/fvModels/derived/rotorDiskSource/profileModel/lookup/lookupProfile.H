/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::lookupProfile

Description
    Look-up based profile data - drag and lift coefficients are linearly
    interpolated based on the supplied angle of attack

    Input in list format:

        data
        (
            (AOA1 Cd1 Cl2)
            (AOA2 Cd2 Cl2)
            ...
            (AOAN CdN CdN)
        );

    where:
        AOA = angle of attack [deg] converted to [rad] internally
        Cd  = drag coefficient
        Cl  = lift coefficient

SourceFiles
    lookupProfile.C

\*---------------------------------------------------------------------------*/

#ifndef lookupProfile_H
#define lookupProfile_H

#include "profileModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class lookupProfile Declaration
\*---------------------------------------------------------------------------*/

class lookupProfile
:
    public profileModel
{

protected:

    // Protected data

        //- List of angle-of-attack values [deg] on input, converted to [rad]
        List<scalar> AOA_;

        //- List of drag coefficient values
        List<scalar> Cd_;

        //- List of lift coefficient values
        List<scalar> Cl_;


    // Protected Member Functions

        //- Return the interpolation indices and gradient
        void interpolateWeights
        (
            const scalar& xIn,
            const List<scalar>& values,
            label& i1,
            label& i2,
            scalar& ddx
        ) const;


public:

    //- Runtime type information
    TypeName("lookup");

    //- Constructor
    lookupProfile(const dictionary& dict, const word& modelName);


    // Member Functions

        // Evaluation

            //- Return the Cd and Cl for a given angle-of-attack
            virtual void Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
