/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionObjects::adjustTimeStepToChemistry

Description
    Returns the minimum chemistry chemical time scale

    This allows the solver to temporally resolve chemical changes, in order to
    better couple the chemistry and transport, or in order to accurately
    post-process the chemical changes.

    Note that this function only does anything if time step adjustment is
    enabled in the controlDict.

    Example of function object specification:
    \verbatim
    adjustTimeStepToChemistry
    {
        type            adjustTimeStepToChemistry;
        libs            ("libchemistryModel.so");
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: adjustTimeStepToChemistry | yes |
        phase        | name of the reacting phase | no       |
    \endtable

SourceFiles
    adjustTimeStepToChemistry.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_adjustTimeStepToChemistry_H
#define functionObjects_adjustTimeStepToChemistry_H

#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class adjustTimeStepToChemistry Declaration
\*---------------------------------------------------------------------------*/

class adjustTimeStepToChemistry
:
    public regionFunctionObject
{
    // Private Data

        //- Name of the phase
        word phaseName_;


public:

    //- Runtime type information
    TypeName("adjustTimeStepToChemistry");


    // Constructors

        //- Construct from components
        adjustTimeStepToChemistry
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        adjustTimeStepToChemistry(const adjustTimeStepToChemistry&) = delete;


    // Destructor
    virtual ~adjustTimeStepToChemistry();


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
        void operator=(const adjustTimeStepToChemistry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
