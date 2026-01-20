/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::radiationModels::absorptionEmission::wideBandReaction

Description
    wideBandReaction radiation absorption and emission coefficients for
    continuous phase. Exactly as wideBand, but with additional effects of the
    reaction heat release rate included.

See also
    Foam::radiationModels::absorptionEmissionModels::wideBand

SourceFiles
    wideBandReaction.C

\*---------------------------------------------------------------------------*/

#ifndef wideBandReaction_H
#define wideBandReaction_H

#include "wideBand.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                     Class wideBandReaction Declaration
\*---------------------------------------------------------------------------*/

class wideBandReaction
:
    public wideBand
{
private:

        //- Proportion of the heat released rate emitted
        FixedList<scalar, maxBands_> iEhrrCoeffs_;


public:

    //- Runtime type information
    TypeName("wideBandReaction");


    // Constructors

        //- Construct from components
        wideBandReaction(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wideBandReaction();


    // Member Functions

        //- Emission contribution for continuous phase
        tmp<volScalarField> ECont(const label bandi = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
