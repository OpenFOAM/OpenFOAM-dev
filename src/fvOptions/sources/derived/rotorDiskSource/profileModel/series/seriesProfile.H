/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::seriesProfile

Description
    Series-up based profile data - drag and lift coefficients computed as
    sum of cosine series

        Cd = sum_i(CdCoeff)*cos(i*AOA)
        Cl = sum_i(ClCoeff)*sin(i*AOA)

    where:
        AOA = angle of attack [deg] converted to [rad] internally
        Cd = drag coefficent
        Cl = lift coefficent

    Input in two (arbitrary length) lists:

        CdCoeffs (coeff1 coeff2 ... coeffN);
        ClCoeffs (coeff1 coeff2 ... coeffN);

SourceFiles
    seriesProfile.C

\*---------------------------------------------------------------------------*/

#ifndef seriesProfile_H
#define seriesProfile_H

#include "profileModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class seriesProfile Declaration
\*---------------------------------------------------------------------------*/

class seriesProfile
:
    public profileModel
{

protected:

    // Protected data

        //- List of drag coefficient values
        List<scalar> CdCoeffs_;

        //- List of lift coefficient values
        List<scalar> ClCoeffs_;


    // Protected Member Functions

        // Evaluate

            //- Drag
            scalar evaluateDrag
            (
                const scalar& xIn,
                const List<scalar>& values
            ) const;

            //- Lift
            scalar evaluateLift
            (
                const scalar& xIn,
                const List<scalar>& values
            ) const;


public:

    //- Runtime type information
    TypeName("series");

    //- Constructor
    seriesProfile(const dictionary& dict, const word& modelName);


    // Member functions

        // Evaluation

            //- Return the Cd and Cl for a given angle-of-attack
            virtual void Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
