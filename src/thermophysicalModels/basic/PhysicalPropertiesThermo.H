/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::PhysicalPropertiesThermo

Description
    Wrapper around a thermo which also constructs the physical properties
    dictionary

\*---------------------------------------------------------------------------*/

#ifndef PhysicalPropertiesThermo_H
#define PhysicalPropertiesThermo_H

#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class PhysicalPropertiesThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermoType>
class PhysicalPropertiesThermo
:
    public physicalProperties,
    public BasicThermoType
{
public:

    // Constructors

        //- Construct from mesh and phase name
        PhysicalPropertiesThermo
        (
            const fvMesh& mesh,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~PhysicalPropertiesThermo();


    // Member Functions

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const;

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties();

        //- Read the properties dictionary
        using BasicThermoType::read;

        //- Read the properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhysicalPropertiesThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
