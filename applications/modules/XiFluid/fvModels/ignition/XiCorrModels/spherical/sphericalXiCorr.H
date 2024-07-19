/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiCorrModels::spherical

Description
    Spherical ignition kernel flame-wrinkling correction model

Usage
    Example usage:
    \verbatim
    XiCorr
    {
        type            spherical;
        select          all;
        sphereFraction  1;
    }
    \endverbatim

    Where:
    \table
        Property       | Description             | Required  | Default value
        select         | Correction cell set selection | yes      |
        sphereFraction | Kernel sphere fraction  | no        | 1
        bMin           | Min b below which no correction | no | 0.01
    \endtable

SourceFiles
    spherical.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalXiCorr_H
#define sphericalXiCorr_H

#include "XiCorrModel.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiCorrModels
{

/*---------------------------------------------------------------------------*\
                           Class spherical Declaration
\*---------------------------------------------------------------------------*/

class spherical
:
    public XiCorrModel
{
    // Private Data

        //- Kernel sphere fraction
        //  If the ignition is adjacent to n symmetry planes: 1/2^n
        dimensionedScalar sphereFraction_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct from mesh and dictionary
        spherical(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        spherical(const spherical&) = delete;


    //- Destructor
    virtual ~spherical();


    // Member Functions

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const spherical&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiCorrModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
