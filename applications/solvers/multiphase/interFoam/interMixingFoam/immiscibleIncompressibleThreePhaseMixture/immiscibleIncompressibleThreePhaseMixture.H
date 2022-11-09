/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::immiscibleIncompressibleThreePhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    immiscibleIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleThreePhaseMixture_H
#define immiscibleIncompressibleThreePhaseMixture_H

#include "incompressibleThreePhaseMixture.H"
#include "threePhaseInterfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class immiscibleIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleThreePhaseMixture
:
    public incompressibleThreePhaseMixture,
    public threePhaseInterfaceProperties
{
public:

    // Constructors

        //- Construct from the mixture velocity
        immiscibleIncompressibleThreePhaseMixture(const volVectorField& U);


    //- Destructor
    virtual ~immiscibleIncompressibleThreePhaseMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            incompressibleThreePhaseMixture::correct();
            threePhaseInterfaceProperties::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
