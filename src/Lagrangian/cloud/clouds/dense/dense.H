/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::clouds::dense

Description
    Base class for clouds with a user-specified density field

SourceFiles
    dense.C

\*---------------------------------------------------------------------------*/

#ifndef dense_H
#define dense_H

#include "massive.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class denseBase Declaration
\*---------------------------------------------------------------------------*/

class denseBase
{
protected:

    // Protected Data

        //- Density
        LagrangianScalarDynamicField rho_;


public:

    // Constructors

        //- Construct from a reference to the cloud
        denseBase(const cloud&);
};

/*---------------------------------------------------------------------------*\
                             Class dense Declaration
\*---------------------------------------------------------------------------*/

class dense
:
    private denseBase,
    public massive
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("dense");


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        dense(const cloud&, const shaped&);


    //- Destructor
    virtual ~dense();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
