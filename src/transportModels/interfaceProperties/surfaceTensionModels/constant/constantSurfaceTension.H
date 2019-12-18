/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::surfaceTensionModels::constant

Description
    Uniform constant surface tension model.

Usage
    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type                constant;
            sigma               0.07;
        }
    \endverbatim

See also
    Foam::surfaceTensionModel

SourceFiles
    constantSurfaceTension.C

\*---------------------------------------------------------------------------*/

#ifndef constantSurfaceTension_H
#define constantSurfaceTension_H

#include "surfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceTensionModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public surfaceTensionModel
{
    // Private Data

        //- Surface tension coefficient
        dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary and mesh
        constant
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Surface tension coefficient
        virtual tmp<volScalarField> sigma() const;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict);

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
