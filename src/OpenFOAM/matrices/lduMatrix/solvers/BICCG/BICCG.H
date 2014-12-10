/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::BICCG

Description
    Diagonal incomplete LU preconditioned BiCG solver derived from the general
    preconditioned BiCG solver PBiCG but with the choice of preconditioner
    pre-selected.

Deprecated
    This solver is present for backward-compatibility and the PBiCG solver
    should be used instead. (deprecated Apr 2008)

SourceFiles
    BICCG.C

\*---------------------------------------------------------------------------*/

#ifndef BICCG_H
#define BICCG_H

#include "PBiCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BICCG Declaration
\*---------------------------------------------------------------------------*/

class BICCG
:
    public PBiCG
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BICCG(const BICCG&);

        //- Disallow default bitwise assignment
        void operator=(const BICCG&);

public:

    //- Return the dictionary constructed from the components.
    //  Needed for backward compatibility
    static dictionary solverDict
    (
        const scalar tol,
        const scalar relTol
    );

    //- Return the dictionary constructed from the old-style data-stream.
    //  Needed for backward compatibility
    static dictionary solverDict(Istream&);


    //- Runtime type information
    TypeName("BICCG");


    // Constructors

        //- Construct from matrix components and solver data stream
        BICCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


        //- Construct from matrix components and tolerances
        BICCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const scalar tolerance,
            const scalar relTol = 0.0
        );

    //- Destructor
    virtual ~BICCG()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
