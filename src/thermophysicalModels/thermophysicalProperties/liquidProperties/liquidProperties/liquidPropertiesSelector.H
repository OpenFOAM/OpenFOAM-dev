/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Foam::liquidPropertiesSelector

Description
    Wrapper class providing run-time selection of thermophysicalProperties
    for the templated thermodynamics packages.

SourceFiles
    liquidPropertiesSelectorI.H
    liquidPropertiesSelector.C

\*---------------------------------------------------------------------------*/

#ifndef liquidPropertiesSelector_H
#define liquidPropertiesSelector_H

#include "thermophysicalPropertiesSelector.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class liquidPropertiesSelector Declaration
\*---------------------------------------------------------------------------*/

class liquidPropertiesSelector
:
    public thermophysicalPropertiesSelector<liquidProperties>
{
public:

    // Constructors

        //- Inherit constructors
        using
            thermophysicalPropertiesSelector<liquidProperties>::
            thermophysicalPropertiesSelector;


    //- Destructor
    virtual ~liquidPropertiesSelector()
    {}


    // Member Functions

        // Physical properties

            //- Surface tension [N/m]
            virtual scalar sigma(const scalar p, const scalar T) const;

            //- Viscosity [Pa s]
            inline scalar mu(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidPropertiesSelectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
