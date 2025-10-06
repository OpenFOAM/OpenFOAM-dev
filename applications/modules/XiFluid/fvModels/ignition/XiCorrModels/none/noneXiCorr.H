/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::XiCorrModels::none

Description
    Dummy ignition kernel flame-wrinkling correction model

    which does not apply any correction

Usage
    Example usage:
    \verbatim
    XiCorr
    {
        type            none;
        cellZone        all;
    }
    \endverbatim
    or simply omit the \c XiCorr entry and the \c none model will be selected
    by default.

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef noneXiCorr_H
#define noneXiCorr_H

#include "XiCorrModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiCorrModels
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public XiCorrModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from mesh and dictionary
        none(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const
        {
            NotImplemented;
            return dimensionedScalar(dimless, scalar(0));
        }

        virtual void XiCorr
        (
            volScalarField& Xi,
            const volScalarField& b,
            const volScalarField& mgb
        ) const
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiCorrModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
