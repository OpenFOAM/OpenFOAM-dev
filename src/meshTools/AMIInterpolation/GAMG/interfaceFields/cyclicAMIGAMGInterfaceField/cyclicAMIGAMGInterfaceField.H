/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::cyclicAMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    cyclicAMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIGAMGInterfaceField_H
#define cyclicAMIGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "cyclicAMIGAMGInterface.H"
#include "cyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicAMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    public cyclicAMILduInterfaceField
{
    // Private Data

        //- Local reference cast into the cyclic interface
        const cyclicAMIGAMGInterface& cyclicAMIInterface_;

        //- Rank of component for transformation
        int rank_;


public:

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        cyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        cyclicAMIGAMGInterfaceField
        (
            const cyclicAMIGAMGInterfaceField&
        ) = delete;


    //- Destructor
    virtual ~cyclicAMIGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicAMIInterface_.size();
            }


        // Interface matrix update

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        //- Cyclic interface functions

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicAMIInterface_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cyclicAMIGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
