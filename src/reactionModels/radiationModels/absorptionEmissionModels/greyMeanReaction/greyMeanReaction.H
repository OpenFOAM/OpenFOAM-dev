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
    Foam::radiationModels::absorptionEmission::greyMeanReaction

Description
    greyMeanReaction radiation absorption and emission coefficients for
    continuous phase. Exactly as greyMean, but with additional effects of the
    reaction heat release rate included.

See also
    Foam::radiationModels::absorptionEmissionModels::greyMean

SourceFiles
    greyMeanReaction.C

\*---------------------------------------------------------------------------*/

#ifndef greyMeanReaction_H
#define greyMeanReaction_H

#include "greyMean.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                     Class greyMeanReaction Declaration
\*---------------------------------------------------------------------------*/

class greyMeanReaction
:
    public greyMean
{
private:

    // Private Data

        //- Emission constant coefficient
        const scalar EhrrCoeff_;


public:

    //- Runtime type information
    TypeName("greyMeanReaction");


    // Constructors

        //- Construct from components
        greyMeanReaction(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~greyMeanReaction();


    // Member Functions

        //- Emission contribution for continuous phase
        tmp<volScalarField> ECont(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
