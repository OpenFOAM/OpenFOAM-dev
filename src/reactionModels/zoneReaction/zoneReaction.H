/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::reactionModels::zoneReaction

Description
    Zone-filtered reaction model.

    Enable the reactions within the specified list of cell-zones and set
    to zero elsewhere.

SourceFiles
    zoneReaction.C

\*---------------------------------------------------------------------------*/

#ifndef zoneReaction_H
#define zoneReaction_H

#include "reactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                       Class zoneReaction Declaration
\*---------------------------------------------------------------------------*/

class zoneReaction
:
    public reactionModel
{
    // Private Data

        //- The reaction model to be zone-filtered
        autoPtr<reactionModel> reactionModelPtr_;

        //- List of zone names in which the reactions are active
        wordList zoneNames_;


    // Private Member Functions

        //- Filter the reaction-rate matrix on the cellZones
        tmp<fvScalarMatrix> filter(const tmp<fvScalarMatrix>& tR) const;

        //- Filter the given field on the cellZones
        template<class GeoField>
        inline tmp<GeoField> filter(const tmp<GeoField>& tS) const;


public:

    //- Runtime type information
    TypeName("zoneReaction");


    // Constructors

        //- Construct from components
        zoneReaction
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& reactionProperties
        );

        //- Disallow default bitwise copy construction
        zoneReaction(const zoneReaction&) = delete;


    //- Destructor
    virtual ~zoneReaction();


    // Member Functions

        //- Correct reaction rate
        virtual void correct();

        //- Specie consumption rate field
        virtual tmp<volScalarField::Internal> R(const label speciei) const;

        //- Specie consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zoneReaction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
