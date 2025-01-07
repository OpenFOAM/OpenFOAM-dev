/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::drag

Description
    Base class for Lagrangian drag models

SourceFiles
    drag.C

\*---------------------------------------------------------------------------*/

#ifndef drag_H
#define drag_H

#include "cloudLagrangianModel.H"
#include "CloudDerivedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class drag Declaration
\*---------------------------------------------------------------------------*/

class drag
:
    public LagrangianModel,
    protected cloudLagrangianModel
{
private:

    // Private Member Functions

        //- Calculate the drag coefficient
        virtual tmp<LagrangianSubScalarField> calcD
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const = 0;

        //- Add a source term to the velocity equation
        void addUSup
        (
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("drag");


    // Public Data

        //- Drag coefficient
        const CloudDerivedField<scalar>& D;


    // Constructors

        //- Construct from components
        drag
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return the name of the velocity field
        virtual wordList addSupFields() const;

        //- Add a source term to the velocity equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;

        //- Add a source term to the momentum equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
