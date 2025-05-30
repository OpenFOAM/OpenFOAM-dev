/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
    Foam::wallBoilingModels::departureDiameterModels::TolubinskiKostanchuk

Description
    Tolubinski-Kostanchuk correlation for bubble departure diameter.

    Reference:
    \verbatim
        Tolubinsky, V. I., & Kostanchuk, D. M. (1970).
        Vapour bubbles growth rate and heat transfer intensity at subcooled
        water boiling.
        In International Heat Transfer Conference 4 (Vol. 23). Begel House Inc.
    \endverbatim

SourceFiles
    TolubinskiKostanchuk.C

\*---------------------------------------------------------------------------*/

#ifndef TolubinskiKostanchuk_H
#define TolubinskiKostanchuk_H

#include "departureDiameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{

/*---------------------------------------------------------------------------*\
                    Class TolubinskiKostanchuk Declaration
\*---------------------------------------------------------------------------*/

class TolubinskiKostanchuk
:
    public departureDiameterModel
{
    // Private Data

        //- Coefficient of the temperature term
        const dimensionedScalar dRef_;

        //- Maximum diameter
        const dimensionedScalar dMax_;

        //- Minimum diameter
        const dimensionedScalar dMin_;


    // Private Member Functions

        //- Calculate and return the bubble departure diameter
        template<class ScalarFieldType>
        tmp<ScalarFieldType> calculate
        (
            const ScalarFieldType& Tl,
            const ScalarFieldType& Tsatw
        ) const;


public:

    //- Runtime type information
    TypeName("TolubinskiKostanchuk");


    // Constructors

        //- Construct from a dictionary
        TolubinskiKostanchuk(const dictionary& dict);

        //- Copy construct
        TolubinskiKostanchuk(const TolubinskiKostanchuk& model);

        //- Construct and return a clone
        virtual autoPtr<departureDiameterModel> clone() const
        {
            return autoPtr<departureDiameterModel>
            (
                new TolubinskiKostanchuk(*this)
            );
        }


    //- Destructor
    virtual ~TolubinskiKostanchuk();


    // Member Functions

        //- Calculate and return the departure diameter field for a wall patch
        virtual tmp<scalarField> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapour,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        //- Calculate and return the departure diameter field for a stationary
        //  phase and fluid interface
        virtual tmp<volScalarField::Internal> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const phaseModel& solid,
            const volScalarField::Internal& Tf,
            const volScalarField::Internal& Tsatw,
            const volScalarField::Internal& L
        ) const;

        //- Calculate and return the departure diameter field for a stationary
        //  phase and fluid interface
        virtual tmp<volScalarField> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapour,
            const phaseModel& solid,
            const volScalarField& Tf,
            const volScalarField& Tsatw,
            const volScalarField& L
        ) const;

        //- Write to stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace departureDiameterModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
