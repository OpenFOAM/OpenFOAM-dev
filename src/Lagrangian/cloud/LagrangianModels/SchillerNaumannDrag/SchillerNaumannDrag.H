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
    Foam::SchillerNaumannDrag

Description
    Drag model of Schiller and Naumann

    References:
    \verbatim
        Schiller, L. & Naumann, A. (1933).
        A drag coefficient correlation.
        Zeit. Ver. Deutsch. Ing., 77, 318-320.
    \endverbatim

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        SchillerNaumannDrag;
    }
    \endverbatim

SourceFiles
    SchillerNaumannDrag.C

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumannDrag_H
#define SchillerNaumannDrag_H

#include "drag.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class SchillerNaumannDrag Declaration
\*---------------------------------------------------------------------------*/

class SchillerNaumannDrag
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


public:

    //- Runtime type information
    TypeName("SchillerNaumannDrag");


    // Constructors

        //- Construct from components
        SchillerNaumannDrag
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return the drag coefficient times Reynold's number, as a function
        //  of the Reynolds number. Used by other models.
        static tmp<LagrangianSubScalarField> CdRe
        (
            const LagrangianSubScalarField& Re
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
