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
    Foam::cyclicGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    cyclicGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGAMGInterfaceField_H
#define cyclicGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "cyclicGAMGInterface.H"
#include "cyclicLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicGAMGInterfaceField
:
    public GAMGInterfaceField,
    public cyclicLduInterfaceField
{
    // Private Data

        //- Local reference cast into the cyclic interface
        const cyclicGAMGInterface& cyclicInterface_;

        //- Rank of component for transformation
        int rank_;


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        cyclicGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const int rank
        );

        //- Disallow default bitwise copy construction
        cyclicGAMGInterfaceField(const cyclicGAMGInterfaceField&) = delete;


    //- Destructor
    virtual ~cyclicGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicInterface_.size();
            }


        //- Cyclic interface functions

            //- Return the transformation
            virtual const transformer& transform() const
            {
                return cyclicInterface_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cyclicGAMGInterfaceField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
