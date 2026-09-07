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
    Foam::fv::sizeSpecificMassTransferIntermediate

Description
    Mix-in interface for size-specific mass-transfer models; i.e., fvModels
    which specify mass-transfer rates for individual size-groups. Provides
    access to the source terms for the groups.

SourceFiles
    sizeSpecificMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef sizeSpecificMassTransfer_H
#define sizeSpecificMassTransfer_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class sizeSpecificMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class sizeSpecificMassTransfer
{
protected:

    // Protected Data

        //- Model name
        const word& name_;

        //- Population balance model
        const populationBalanceModel& popBal_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the total phase change rate
        tmp<volInternalScalarField> mDot() const;


public:

    //- Runtime type information
    TypeName("sizeSpecificMassTransfer");


    // Constructors

        //- Construct from a mesh and coefficients
        sizeSpecificMassTransfer
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor
    virtual ~sizeSpecificMassTransfer();


    // Member Functions

        //- Get the group index given the field name, or else -1
        label groupi(const word& fieldName) const;

        //- Return the phase change rate divided by the group fraction for
        //  a given group index
        virtual tmp<volInternalScalarField> groupMDotByF(const label) const = 0;

        //- Return the phase change rate for a given group
        virtual tmp<volInternalScalarField> groupMDot(const label) const;

        //- Return the source for a given group index
        virtual tmp<volInternalScalarField> groupS(const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
