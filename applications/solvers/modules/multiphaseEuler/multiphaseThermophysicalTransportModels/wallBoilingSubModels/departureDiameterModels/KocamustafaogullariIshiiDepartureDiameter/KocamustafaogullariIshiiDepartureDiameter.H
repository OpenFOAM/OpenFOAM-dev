/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::wallBoilingModels::departureDiameterModels::
    KocamustafaogullariIshiiDepartureDiameter

Description
    A correlation for bubble departure diameter.

    Requires model parameter 'phi': contact angle in degrees.

    Reference:
    \verbatim
        Kocamustafaogullari, G., & Ishii, M. (1983).
        Interfacial area and nucleation site density in boiling systems.
        International Journal of Heat and Mass Transfer, 26(9), 1377-1387.
    \endverbatim

SourceFiles
    KocamustafaogullariIshiiDepartureDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef KocamustafaogullariIshiiDepartureDiameter_H
#define KocamustafaogullariIshiiDepartureDiameter_H

#include "departureDiameterModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{

/*---------------------------------------------------------------------------*\
          Class KocamustafaogullariIshiiDepartureDiameter Declaration
\*---------------------------------------------------------------------------*/

class KocamustafaogullariIshiiDepartureDiameter
:
    public departureDiameterModel
{
    // Private Data

        //- Contact angle
        const dimensionedScalar phi_;


    // Private Member Functions

        //- Calculate and return the bubble departure diameter
        template<class ScalarFieldType>
        tmp<ScalarFieldType> calculate
        (
            const fvMesh& mesh,
            const ScalarFieldType& Tl,
            const ScalarFieldType& Tsatw,
            const ScalarFieldType& L,
            const ScalarFieldType& rhoLiquid,
            const ScalarFieldType& rhoVapor,
            const ScalarFieldType& sigma
        ) const;


public:

    //- Runtime type information
    TypeName("KocamustafaogullariIshii");


    // Constructors

        //- Construct from a dictionary
        KocamustafaogullariIshiiDepartureDiameter(const dictionary& dict);

        //- Copy construct
        KocamustafaogullariIshiiDepartureDiameter
        (
            const KocamustafaogullariIshiiDepartureDiameter& model
        );

        //- Construct and return a clone
        virtual autoPtr<departureDiameterModel> clone() const
        {
            return autoPtr<departureDiameterModel>
            (
                new KocamustafaogullariIshiiDepartureDiameter(*this)
            );
        }


    //- Destructor
    virtual ~KocamustafaogullariIshiiDepartureDiameter();


    // Member Functions

        //- Calculate and return the departure diameter field for a wall patch
        virtual tmp<scalarField> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        //- Calculate and return the departure diameter field for a stationary
        //  phase and fluid interface
        virtual tmp<volScalarField> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
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
