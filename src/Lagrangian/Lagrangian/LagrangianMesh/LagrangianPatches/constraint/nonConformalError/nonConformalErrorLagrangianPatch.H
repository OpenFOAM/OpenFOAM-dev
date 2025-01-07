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
    Foam::nonConformalErrorLagrangianPatch

Description
    Non-conformal-error Lagrangian patch

SourceFiles
    nonConformalErrorLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorLagrangianPatch_H
#define nonConformalErrorLagrangianPatch_H

#include "LagrangianPatch.H"
#include "nonConformalErrorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nonConformalErrorLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalErrorLagrangianPatch
:
    public LagrangianPatch
{
private:

    // Private Data

        //- Reference to the nonConformalError poly patch
        const nonConformalErrorPolyPatch& nonConformalErrorPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalErrorPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        nonConformalErrorLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~nonConformalErrorLagrangianPatch();


    // Member Functions

        //- Return the nonConformalError patch
        const nonConformalErrorPolyPatch& nonConformalErrorPatch() const
        {
            return nonConformalErrorPatch_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
