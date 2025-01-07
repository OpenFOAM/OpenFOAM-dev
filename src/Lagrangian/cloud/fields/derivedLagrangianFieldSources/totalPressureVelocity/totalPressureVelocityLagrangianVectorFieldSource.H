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
    Foam::totalPressureVelocityLagrangianVectorFieldSource

Description
    This source condition provides a velocity in a given direction with a
    magnitude calculated from the difference between a specified upstream
    pressure and the local interpolated pressure.

Usage
    \table
        Property  | Description                           | Required? | Default
        direction | The velocity direction                | yes       |
        p0        | The upstream total pressure           | yes       |
        pc        | The name of the carrier pressure field | no       | p
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            totalPressureVelocity;
        direction       (1 0 0);
        p0              1 [bar];
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    totalPressureVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureVelocityLagrangianVectorFieldSource_H
#define totalPressureVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "CloudLagrangianFieldSource.H"
#include "Function1LagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class totalPressureVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class totalPressureVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private Function1LagrangianFieldSource<vector>,
    private CloudLagrangianFieldSource<vector>
{
private:

    // Private Data

        //- Cached dictionary. This is needed as construction of the total
        //  pressure function is delayed until the dimensions of pressure are
        //  known. This will be thrown away once total pressure has been
        //  constructed.
        mutable autoPtr<dictionary> dict_;

        //- Direction
        const autoPtr<Function1<vector>> direction_;

        //- Total pressure
        mutable autoPtr<Function1<scalar>> p0_;

        //- Name of the carrier density field
        const word rhocName_;

        //- Name of the carrier pressure field
        const word pcName_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource<vector>::value;


public:

    //- Runtime type information
    TypeName("totalPressureVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        totalPressureVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        totalPressureVelocityLagrangianVectorFieldSource
        (
            const totalPressureVelocityLagrangianVectorFieldSource&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianVectorFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianVectorFieldSource>
            (
                new totalPressureVelocityLagrangianVectorFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~totalPressureVelocityLagrangianVectorFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubVectorField> value
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
