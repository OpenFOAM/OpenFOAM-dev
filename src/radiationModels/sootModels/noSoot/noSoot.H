/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::radiationModels::sootModels::noSoot

Description
    Dummy soot model for 'none'

SourceFiles
    noSoot.C

\*---------------------------------------------------------------------------*/

#ifndef noSoot_H
#define noSoot_H

#include "sootModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{

/*---------------------------------------------------------------------------*\
                           Class noSoot Declaration
\*---------------------------------------------------------------------------*/

class noSoot
:
    public sootModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noSoot(const dictionary& dict, const fvMesh& mesh, const word&);


    //- Destructor
    virtual ~noSoot();


    // Member Functions

        // Edit

            //- Main update/correction routine
            void correct();

        // Access

            //- Return Ysoot
            const volScalarField& soot() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sootModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
