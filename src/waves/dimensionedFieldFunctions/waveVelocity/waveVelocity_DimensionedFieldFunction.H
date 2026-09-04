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
    Foam::DimensionedFieldFunctions::waveVelocity

Description
    Initialises the velocity to that specified by a superposition of wave
    models. All the wave modelling parameters are obtained from a centrally
    registered waveSuperposition class.

Usage
    Example of the boundary condition specification:
    \verbatim
    internalField
    {
        type        waveVelocity;
        libs        ("libwaves.so");
    }
    \endverbatim

See also
    Foam::waveVelocityFvPatchScalarField
    Foam::waveSuperposition

SourceFiles
    waveVelocity_DimensionedFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef waveVelocity_DimensionedFieldFunction_H
#define waveVelocity_DimensionedFieldFunction_H

#include "DimensionedFieldFunction.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace DimensionedFieldFunctions
{

/*---------------------------------------------------------------------------*\
                            Class waveVelocity Declaration
\*---------------------------------------------------------------------------*/

class waveVelocity
:
    public DimensionedFieldFunction<DimensionedField<vector, fvMesh>>
{

public:

    //- Runtime type information
    TypeName("waveVelocity");


    // Constructors

        //- Construct with dictionary to initialise given field
        waveVelocity
        (
            const dictionary& dict,
            DimensionedField<vector, fvMesh>& field
        );

        //- Construct a copy for the given field
        waveVelocity
        (
            const waveVelocity& dff,
            DimensionedField<vector, fvMesh>& field
        );

        //- Construct and return a clone for the specified field
        virtual autoPtr
        <
            DimensionedFieldFunction<DimensionedField<vector, fvMesh>>
        > clone(DimensionedField<vector, fvMesh>& field) const;


    //- Destructor
    virtual ~waveVelocity()
    {}


    // Member Functions

        //- Evaluate the function and set the field
        virtual void evaluate();

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DimensionedFieldFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
