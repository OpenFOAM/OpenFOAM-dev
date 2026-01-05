/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::Lagrangian::heatTransfer

Description
    Base class for Lagrangian heat transfer models

SourceFiles
    heatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransfer_H
#define heatTransfer_H

#include "cloudLagrangianModel.H"
#include "thermal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class heatTransfer Declaration
\*---------------------------------------------------------------------------*/

class heatTransfer
:
    public LagrangianModel,
    protected cloudLagrangianModel
{
private:

    // Private Member Data

        //- Thermal cloud reference
        const clouds::thermal& thermalCloud_;


    // Private Member Functions

        //- Calculate the heat transfer coefficient
        virtual tmp<LagrangianSubScalarField> calcH
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const = 0;


public:

    //- Runtime type information
    TypeName("heatTransfer");


    // Public Data

        //- Heat transfer coefficient
        const CloudDerivedField<scalar>& H;


    // Constructors

        //- Construct from components
        heatTransfer
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return the name of the energy field
        virtual wordList addSupFields() const;

        //- Return true for the energy or carrier energy field
        virtual bool addsSupToField
        (
            const word& fieldName,
            const word& eqnFieldName
        ) const;

        //- Add a source term to the energy equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubScalarSubField& e,
            LagrangianEqn<scalar>& eqn
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
