/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::functionObjects::log

Description
    Calculates the natural logarithm of the specified scalar field.

    Performs \f$ln(max(x, a))\f$ where \f$x\f$ is the field and \f$a\f$ an
    optional clip to handle 0 or negative \f$x\f$.  Dimension checking can
    optionally be suspended for this operation if \f$x\f$ is dimensioned.

    Example of function object specification:
    \verbatim
    log1
    {
        type            log;
        libs            ("libfieldFunctionObjects.so");
        field           p;
        clip            1e-3;
        checkDimensions no;
    }
    \endverbatim

    or using \c foamPostProcess
    \verbatim
    foamPostProcess -func 'log(p, clip=1e-3, checkDimensions=no)'
    \endverbatim

Usage
    \table
        Property        | Description               | Required | Default value
        type            | Type name: log            | yes      |
        clip            | Clip value                | no       |
        checkDimensions | Dimension checking switch | no       |
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    log.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_log_H
#define functionObjects_log_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                             Class log Declaration
\*---------------------------------------------------------------------------*/

class log
:
    public fieldExpression
{
    // Private Data

        //- Optional clip, default = false
        bool clip_;

        //- Optional clip value
        scalar clipValue_;

        //- Optional dimension checking override
        Switch checkDimensions_;


    // Private Member Functions

        //- Calculate the log field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("log");


    // Constructors

        //- Construct from Time and dictionary
        log
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~log();


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
