/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::power

Description
    Evaluates and writes the following power fields.
    \verbatim
    + tauUSf:           tauf & Uf
    + divTauU:          fvc::div(tauUSf)
    + divTauDotU:       fvc::div(tauf) & U
    + tauDDotGradU:     divTauU - divTauDotU

    + pUSf:            -pIf & Uf
    + divPU:            fvc::div(pUSf)
    + gradPU:          -fvc::grad(p) & U
    + pDivU:            divPU - gradPU

    + stressUSf:        tauUSf + pUSf
    + divStressU:       fvc::div(stressUSf)
    + divStressDotU:    fvc::div(stressf) & U
    + stressDDotGradU:  divStressU - divStressDotU
    \endverbatim

    where:
    \verbatim
    + p, U = cell centre pressure and velocity
    + pIf, Uf = face normal pressure force, and velocity, vectors
    + tauf = face shear force vectors [-transport.devTau()*mesh().magSf()]
    + stressf = tauf - pIf
    \endverbatim

    The domain integral of the volume fields (all except '...Sf' fields) is also
    printed to a log file in the postProcessing directory.  Domain integrals are
    scaled by a factor, which defaults to 1e9, to make them easier to process.

Usage
    An example is given below of a function object specification in the
    system/functions file with a scaling factor of 1e6 for the domain integrals
    \verbatim
    power
    {
        type        power;
        libs        ("libfieldFunctionObjects.so");
        executeControl  writeTime;
        writeControl    writeTime;
        log             true;
        factor          1e6;
    }
    \endverbatim

    Alternatively it can be specified using the standard 'power' function object
    configuration file as follows
    \verbatim
    #includeFunc power(factor=1e6)
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    power.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_power_H
#define functionObjects_power_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class power Declaration
\*---------------------------------------------------------------------------*/

class power
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Static data

        //- List of the power field names
        static const wordList fields_;


    // Private Member data

        //- Scaling factor for the power output (defaults to 1e9)
        scalar factor_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("power");


    // Constructors

        //- Construct from Time and dictionary
        power
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        power(const power&) = delete;


    //- Destructor
    virtual ~power();


    // Member Functions

        //- Read the power data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the power field
        virtual bool execute();

        //- Write the power field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const power&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
