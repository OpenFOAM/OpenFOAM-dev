/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::radiationModels::absorptionEmissionModels::cloud

Description
    Retrieves absorption/emission data from a cloud object

SourceFiles
    cloudAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef cloudAbsorptionEmission_H
#define cloudAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public absorptionEmissionModel
{
    // Private Data

        //- Cloud name(s)
        const wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloud");


    // Constructors

        //- Construct from components
        cloud(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~cloud();


    // Member Functions

        //- Absorption coefficient for dispersed phase
        virtual tmp<volScalarField> aDisp(const label bandI = 0) const;

        //- Emission coefficient for dispersed phase
        virtual tmp<volScalarField> eDisp(const label bandI = 0) const;

        //- Return emission contribution for dispersed phase
        virtual tmp<volScalarField> EDisp(const label bandI = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace absorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
