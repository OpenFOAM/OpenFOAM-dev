/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::regionCoupledGAMGInterfaceField

Description
    GAMG agglomerated region coupled interface field.

SourceFiles
    regionCoupledGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledGAMGInterfaceField_H
#define regionCoupledGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "regionCoupledGAMGInterface.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class regionCoupledGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledGAMGInterfaceField
:
    public GAMGInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const regionCoupledGAMGInterface& regionCoupledGAMGInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionCoupledGAMGInterfaceField
        (
            const regionCoupledGAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const regionCoupledGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("regionCoupled");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        regionCoupledGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        regionCoupledGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );

    //- Destructor
    virtual ~regionCoupledGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return regionCoupledGAMGInterface_.size();
            }


            //- Interface matrix update
            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
