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
    Foam::functionObjects::adjustTimeStepToNucleation

Description
    Controls the time-step to the minimum nucleation time scale. A specific
    nucleation model can be specified, otherwise it will scan for and use all
    the nucleation fvModels in operation.

    Note that this function only does anything if time step adjustment is
    enabled in the controlDict.

    Example of function object specification:
    \verbatim
    adjustTimeStepToNucleation
    {
        type            adjustTimeStepToNucleation;
        libs            ("libmultiphaseEulerFunctionObjects.so");
        maxCo           0.1;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: adjustTimeStepToNucleation | yes |
        model        | name of the nucleation fvModel | no   |
        maxCo        | maximum nucleation "Courant" number | no | 1
    \endtable

SourceFiles
    adjustTimeStepToNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_adjustTimeStepToNucleation_H
#define functionObjects_adjustTimeStepToNucleation_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class adjustTimeStepToNucleation Declaration
\*---------------------------------------------------------------------------*/

class adjustTimeStepToNucleation
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the nucleation fvModel
        word modelName_;

        //- Maximum nucleation "Courant" number
        scalar maxCo_;


public:

    //- Runtime type information
    TypeName("adjustTimeStepToNucleation");


    // Constructors

        //- Construct from components
        adjustTimeStepToNucleation
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        adjustTimeStepToNucleation(const adjustTimeStepToNucleation&) = delete;


    // Destructor
    virtual ~adjustTimeStepToNucleation();


    // Member Functions

        //-  Read and reset the timeStep Function1
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Do nothing
        virtual bool execute();

        //- Do nothing
        virtual bool write();

        //- Return the minimum chemistry.deltaTChem()
        virtual scalar maxDeltaT() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adjustTimeStepToNucleation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
