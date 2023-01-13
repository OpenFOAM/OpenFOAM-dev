/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::compressible::alphatPhaseChangeWallFunctionBase

Description
    Abstract base-class for all alphatWallFunctions supporting phase-change.

SourceFiles
    alphatPhaseChangeWallFunctionBase.C

\*---------------------------------------------------------------------------*/

#ifndef alphatPhaseChangeWallFunctionBase_H
#define alphatPhaseChangeWallFunctionBase_H

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseInterface;

namespace compressible
{

/*---------------------------------------------------------------------------*\
             Class alphatPhaseChangeWallFunctionBase Declaration
\*---------------------------------------------------------------------------*/

class alphatPhaseChangeWallFunctionBase
{
protected:

    // Protected data

        //- Name of the phase
        const word phaseName_;

        //- Name of the other phase
        const word otherPhaseName_;


public:

    //- Runtime type information
    TypeName("compressible::alphatPhaseChangeWallFunctionBase");


    // Constructors

        //- Construct null
        alphatPhaseChangeWallFunctionBase();

        //- Construct from a patch, an internal field and a dictionary
        alphatPhaseChangeWallFunctionBase
        (
            const fvPatch& p,
            const DimensionedField<scalar, volMesh>& iF,
            const dictionary&
        );


    // Destructor
    virtual ~alphatPhaseChangeWallFunctionBase();


    // Member Functions

        //- Is there phase change mass transfer for this interface?
        bool activeInterface(const phaseInterface&) const;

        //- Return the rate of phase-change
        virtual const scalarField& dmdtf() const = 0;

        //- Return the rate of phase-change for the given interface
        const scalarField& dmdtf(const phaseInterface&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
