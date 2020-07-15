/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::dragModels::SyamlalOBrien

Description

    Reference:
    \verbatim
        Syamlal, M., Rogers, W., & O’Brien, T. J. (1993).
        Mfix documentation: Theory guide, technical note, doe/metc-94/1004,
        ntis/de94000087.
        National Technical Information Service.
    \endverbatim

SourceFiles
    SyamlalOBrien.C

\*---------------------------------------------------------------------------*/

#ifndef SyamlalOBrien_H
#define SyamlalOBrien_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class SyamlalOBrien Declaration
\*---------------------------------------------------------------------------*/

class SyamlalOBrien
:
    public dragModel
{
public:

    //- Runtime type information
    TypeName("SyamlalOBrien");


    // Constructors

        //- Construct from a dictionary and a phase pair
        SyamlalOBrien
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~SyamlalOBrien();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
