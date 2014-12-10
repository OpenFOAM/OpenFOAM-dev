/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::incompressible::DESModel

Group
    grpIcoDESTurbulence

Description
    Extension of LES models to provide an interface for Detached Eddy
    Simulation turbulence models

SeeAlso
    Foam::LESModel

SourceFiles
    DESModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleDESModel_H
#define incompressibleDESModel_H

// note: use full(er) path so that derived utilities can bring in
// incompressible and compressible variants

#include "incompressible/LES/LESModel/LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                          Class DESModel Declaration
\*---------------------------------------------------------------------------*/

class DESModel
:
    public LESModel
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DESModel(const DESModel&);

        //- Disallow default bitwise assignment
        DESModel& operator=(const DESModel&);


public:

    // Constructors

        //- Construct from components
        DESModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName
        );


    //- Destructor
    virtual ~DESModel()
    {}


    // Public Member Functions

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
