/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::fv::fixedTemperatureConstraint

Description
    Fixed temperature equation constraint

Usage
    \verbatim
    fixedTemperature
    {
        type            fixedTemperatureConstraint;

        selectionMode   all;

        phase           gas; // Optional phase name

        // Uniform temperature constraint
        mode            uniform;
        temperature     constant 500; // Uniform temperature

        // // Looked-up field temperature constraint
        // T            T; // Temperature field name
    }
    \endverbatim

Note
    The 'uniform' option allows the use of a time-varying uniform temperature
    by means of the Function1 type.

SourceFiles
    fixedTemperatureConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureConstraint_H
#define fixedTemperatureConstraint_H

#include "fvConstraint.H"
#include "fvCellSet.H"
#include "NamedEnum.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class fixedTemperatureConstraint Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureConstraint
:
    public fvConstraint
{
public:

    //- Temperature mode
    enum class temperatureMode
    {
        uniform,
        lookup
    };

    //- String representation of mode enums
    static const NamedEnum<temperatureMode, 2> modeNames_;


private:

    // Private Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Operation mode
        temperatureMode mode_;

        //- Uniform temperature [K]
        autoPtr<Function1<scalar>> TValue_;

        //- Temperature field name
        word TName_;

        //- Optional phase name
        word phaseName_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("fixedTemperatureConstraint");


    // Constructors

        //- Construct from components
        fixedTemperatureConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        fixedTemperatureConstraint(const fixedTemperatureConstraint&) = delete;


    //- Destructor
    virtual ~fixedTemperatureConstraint()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Constrain energy equation to fix the temperature
        virtual bool constrain
        (
            fvMatrix<scalar>& eqn,
            const word& fieldName
        ) const;

        //- Update for mesh changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fixedTemperatureConstraint&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
