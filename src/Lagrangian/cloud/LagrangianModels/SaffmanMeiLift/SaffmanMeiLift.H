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
    Foam::SaffmanMeiLift

Description
    Lift model of Saffman (1965) as extended by Mei (1992). Applicable to
    spherical particles.

    References:
    \verbatim
        Saffman, P. G. T. (1965).
        The lift on a small sphere in a slow shear flow.
        Journal of fluid mechanics, 22(2), 385-400.

        Mei, R. (1992).
        An approximate expression for the shear lift force on a spherical
        particle at finite Reynolds number.
        International Journal of Multiphase Flow, 18(1), 145-147.
    \endverbatim

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        SaffmanMeiLift;
    }
    \endverbatim

SourceFiles
    SaffmanMeiLift.C

\*---------------------------------------------------------------------------*/

#ifndef SaffmanMeiLift_H
#define SaffmanMeiLift_H

#include "cloudLagrangianModel.H"
#include "CloudDerivedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class SaffmanMeiLift Declaration
\*---------------------------------------------------------------------------*/

class SaffmanMeiLift
:
    public LagrangianModel,
    private cloudLagrangianModel
{
private:

    // Private Member Functions

        //- Calculate the lift coefficient
        tmp<LagrangianSubTensorField> calcL
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Add a source term to the velocity equation
        void addUSup
        (
            const LagrangianSubVectorSubField& U,
            LagrangianEqn<vector>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("SaffmanMeiLift");


    // Public Data

        //- Lift force
        const CloudDerivedField<tensor>& L;


    // Constructors

        //- Construct from components
        SaffmanMeiLift
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
