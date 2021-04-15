/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef compileThermo_H
#define compileThermo_H

#include "CodedBase.H"
#include "dynamicCode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class compileThermo Declaration
\*---------------------------------------------------------------------------*/

class compileThermo
:
    public CodedBase<compileThermo>
{
    // Private Member Data

        //- Name of the class Thermo basicThermo is instantiated on
        const word thermoTypeName_;

        //- Reference to the dictionary containing the thermo components
        const dictionary& thermoTypeDict_;


    // Private Member Functions

        dictionary optionsDict
        (
            const word& thermoTypeName,
            const word& instantiatedThermoName
        ) const;

        void setFilterVariable
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context,
            const word& name
        ) const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const
        {}


public:

    // Constructors

        //- Construct from name and dictionary
        compileThermo
        (
            const word& thermoTypeName,
            const word& instantiatedThermoName,
            const dictionary& thermoTypeDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
