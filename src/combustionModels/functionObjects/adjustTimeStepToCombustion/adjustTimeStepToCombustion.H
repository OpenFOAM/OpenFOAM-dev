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
    Foam::functionObjects::adjustTimeStepToCombustion

Description
    Adjusts the time step to match bulk reaction time scales. This allows the
    solver to temporally resolve chemical changes, in order to better couple
    the chemistry and transport, or to improve the time-accuracy of
    post-processing.

    Note that this function only does anything if time step adjustment is
    enabled in the controlDict.

    Example of function object specification:
    \verbatim
    adjustTimeStepToCombustion1
    {
        type            adjustTimeStepToCombustion;
        libs            ("libcombustionModels.so");
        maxCo           0.1;
        extrapolate     no;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: adjustTimeStepToCombustion | yes |
        phase        | name of the reacting phase | no       |
        maxCo        | maximum combustion "Courant" number | no | 1
        extrapolate  | try to extrapolate decreases | no     | false
    \endtable

SourceFiles
    adjustTimeStepToCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_adjustTimeStepToCombustion_H
#define functionObjects_adjustTimeStepToCombustion_H

#include "regionFunctionObject.H"
#include "timeIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class adjustTimeStepToCombustion Declaration
\*---------------------------------------------------------------------------*/

class adjustTimeStepToCombustion
:
    public regionFunctionObject
{
    // Private Data

        //- Name of the phase
        word phaseName_;

        //- Maximum combustion "Courant" number. The maximum fraction of the
        //  total combustion time-scale that should be evolved in a single
        //  time-step.
        scalar maxCo_;

        //- Extrapolate reductions in time-step to try and better catch the
        //  "onset" of reaction. Can lead to time-step oscillation.
        bool extrapolate_;

        //- Do we have the previous combustion time-step?
        bool haveCombustionDeltaT0_;

        //- The previous combustion time-step
        scalar combustionDeltaT0_;


    // Private Member Functions

        //- Get the properties dictionary IO object
        typeIOobject<timeIOdictionary> propsDictIo
        (
            const IOobject::readOption& r
        ) const;

public:

    //- Runtime type information
    TypeName("adjustTimeStepToCombustion");


    // Constructors

        //- Construct from components
        adjustTimeStepToCombustion
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        adjustTimeStepToCombustion(const adjustTimeStepToCombustion&) = delete;


    // Destructor
    virtual ~adjustTimeStepToCombustion();


    // Member Functions

        //-  Read and reset the timeStep Function1
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Reset the timeStep from the Function1 of time
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const adjustTimeStepToCombustion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
