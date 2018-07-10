/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::radiation::blackBodyEmission

Description
    Class black body emission

    Table of black body emissive power from:
        Modest, "Radiative Heat Transfer", pp.775-777, 1993

SourceFiles
    blackBodyEmission.C

\*---------------------------------------------------------------------------*/

#ifndef blackBodyEmission_H
#define blackBodyEmission_H

#include "volFields.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                      Class blackBodyEmission Declaration
\*---------------------------------------------------------------------------*/

class blackBodyEmission
{
public:

    //- Static table of black body emissive power
    static const List<Tuple2<scalar, scalar>> emissivePowerTable;


private:

    // Private data

        //- Interpolation table of black body emissive power
        mutable interpolationTable<scalar> table_;

        //- Constant C1
        const dimensionedScalar C1_;

        //- Constant C2
        const dimensionedScalar C2_;

        // Ptr List for black body emission energy field for each wavelength
        PtrList<volScalarField> bLambda_;

        // Reference to the temperature field
        const volScalarField& T_;


    // Private Member Functions

        scalar fLambdaT(const scalar lambdaT) const;


public:

    // Constructors

        //- Construct from components
        blackBodyEmission
        (
            const label nLambda,
            const volScalarField& T
        );


    //- Destructor
    ~blackBodyEmission();


    // Member functions

        // Access

            //- Black body spectrum
            inline const volScalarField& bLambda(const label lambdaI) const
            {
                return bLambda_[lambdaI];
            }

            //- Spectral emission for the black body at T and lambda
            inline dimensionedScalar EblambdaT
            (
                const dimensionedScalar& T,
                const scalar lambda
            ) const
            {
                return (C1_/(pow5(lambda)*(exp(C2_/(lambda*T)) - 1.0)));
            }

            //- Proportion of total energy at T from lambda1 to lambda2
            tmp<Foam::volScalarField> deltaLambdaT
            (
                const volScalarField& T,
                const Vector2D<scalar>& band
            ) const;

            //- Integral energy at T from lambda1 to lambda2
            tmp<Foam::volScalarField> EbDeltaLambdaT
            (
                const volScalarField& T,
                const Vector2D<scalar>& band
            ) const;

    // Edit

        // Update black body emission
        void correct(const label lambdaI, const Vector2D<scalar>& band);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
