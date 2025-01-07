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
    Foam::GidaspowErgunWenYuDrag

Description
    Drag model of Gidaspow which switches between the drag models of Ergun and
    Wen and Yu at different particulate volume dractions

    References:
    \verbatim
        Gidaspow, D. (1994).
        Multiphase flow and fluidisation: continuum and kinetic theory
        descriptions.
        Academic press, New York.
    \endverbatim

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        GidaspowErgunWenYuDrag;
    }
    \endverbatim

SourceFiles
    GidaspowErgunWenYuDrag.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowErgunWenYuDrag_H
#define GidaspowErgunWenYuDrag_H

#include "drag.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class GidaspowErgunWenYuDrag Declaration
\*---------------------------------------------------------------------------*/

class GidaspowErgunWenYuDrag
:
    public drag
{
private:

    // Private Member Functions

        //- Calculate the drag coefficient
        virtual tmp<LagrangianSubScalarField> calcD
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;


    // Private Data

        //- The upper bound on the volume fraction
        const scalar alphaMax_;


public:

    //- Runtime type information
    TypeName("GidaspowErgunWenYuDrag");


    // Constructors

        //- Construct from components
        GidaspowErgunWenYuDrag
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
