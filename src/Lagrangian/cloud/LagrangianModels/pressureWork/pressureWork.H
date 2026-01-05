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
    Foam::Lagrangian::pressureWork

Description
    Model for pressure work. Subtracts the pressure-volume work associated with
    particle expansion from the particle energy, and adds it to the carrier
    energy. This model is needed to make (e.g.,) bubbles reduce in temperature
    as they expand.

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        pressureWork;
    }
    \endverbatim

SourceFiles
    pressureWork.C

\*---------------------------------------------------------------------------*/

#ifndef pressureWork_H
#define pressureWork_H

#include "cloudLagrangianModel.H"
#include "thermal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class pressureWork Declaration
\*---------------------------------------------------------------------------*/

class pressureWork
:
    public LagrangianModel,
    private cloudLagrangianModel
{
private:

    // Private Member Data

        //- Thermal cloud reference
        const clouds::thermal& thermalCloud_;

        //- Pressure at the previous iteration
        autoPtr<LagrangianSubScalarField> pPrevPtr_;


public:

    //- Runtime type information
    TypeName("pressureWork");


    // Constructors

        //- Construct from components
        pressureWork
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

        //- Hook before source evaluation
        virtual void preAddSup
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Add a source term to the energy equation
        virtual void addSup
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& vOrM,
            const LagrangianSubScalarSubField& e,
            LagrangianEqn<scalar>& eqn
        ) const;

        //- Hook after source evaluation
        virtual void postAddSup
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
