/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void kinematicSingleLayer::constrainFilmField
(
    Type& field,
    const typename Type::cmptType& value
)
{
    forAll(intCoupledPatchIDs_, i)
    {
        label patchI = intCoupledPatchIDs_[i];
        field.boundaryField()[patchI] = value;
        if (debug)
        {
            Info<< "Constraining " << field.name()
                << " boundary " << field.boundaryField()[patchI].patch().name()
                << " to " << value << endl;
        }
    }
    forAll(passivePatchIDs_, i)
    {
        label patchI = passivePatchIDs_[i];
        field.boundaryField()[patchI] = value;
        if (debug)
        {
            Info<< "Constraining " << field.name()
                << " boundary " << field.boundaryField()[patchI].patch().name()
                << " to " << value << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace Foam
} // end namespace regionModels
} // end namespace surfaceFilmModels

// ************************************************************************* //
