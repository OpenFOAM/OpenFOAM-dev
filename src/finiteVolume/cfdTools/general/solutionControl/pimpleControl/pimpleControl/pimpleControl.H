/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::pimpleControl

Description

SourceFiles
    pimpleControl.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleControl_H
#define pimpleControl_H

#include "pimpleSingleRegionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pimpleControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleControl
:
    public pimpleNoLoopControl,
    public pimpleSingleRegionControl
{
public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimpleControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        pimpleControl(fvMesh& mesh, const word& algorithmName="PIMPLE");


    //- Destructor
    virtual ~pimpleControl();


    // Member Functions

        using pimpleSingleRegionControl::loop;

        using pimpleSingleRegionControl::run;

        //- Flag to indicate the first iteration
        bool firstPimpleIter() const
        {
            return pimpleSingleRegionControl::firstIter();
        }

        //- Flag to indicate the final iteration
        bool finalPimpleIter() const
        {
            return pimpleSingleRegionControl::finalIter();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
