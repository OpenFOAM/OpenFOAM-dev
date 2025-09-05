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
    Foam::coneDiskDirectionLagrangianVectorFieldSource

Description
    Mix-in for source conditions that provides a random injection direction
    on a disk between two cone angles

SourceFiles
    coneDiskDirectionLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef coneDiskDirectionLagrangianVectorFieldSource_H
#define coneDiskDirectionLagrangianVectorFieldSource_H

#include "LagrangianInjection.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianFieldSourceBase;

/*---------------------------------------------------------------------------*\
            Class coneDiskDirectionLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class coneDiskDirectionLagrangianVectorFieldSource
{
private:

    // Private Data

        //- Reference to the field source
        const LagrangianFieldSourceBase& field_;

        //- Inner half-cone angle
        const autoPtr<Function1<scalar>> thetaInner_;

        //- Outer half-cone angle
        const autoPtr<Function1<scalar>> thetaOuter_;


public:

    // Constructors

        //- Construct from a dictionary
        coneDiskDirectionLagrangianVectorFieldSource
        (
            const LagrangianFieldSourceBase&,
            const dictionary& dict
        );

        //- Copy constructor
        coneDiskDirectionLagrangianVectorFieldSource
        (
            const coneDiskDirectionLagrangianVectorFieldSource&,
            const LagrangianFieldSourceBase&
        );


    //- Destructor
    ~coneDiskDirectionLagrangianVectorFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        tmp<LagrangianSubVectorField> direction
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
