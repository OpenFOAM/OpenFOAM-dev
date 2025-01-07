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
    Foam::symmetryPlaneLagrangianPatch

Description
    Symmetry-plane Lagrangian patch. This is used for representing
    reflected/mirror-image repetitions of the geometry.

SourceFiles
    symmetryPlaneLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneLagrangianPatch_H
#define symmetryPlaneLagrangianPatch_H

#include "LagrangianPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class symmetryPlaneLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class symmetryPlaneLagrangianPatch
:
    public LagrangianPatch
{
public:

    //- Runtime type information
    TypeName(symmetryPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        symmetryPlaneLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~symmetryPlaneLagrangianPatch();


    // Member Functions

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
