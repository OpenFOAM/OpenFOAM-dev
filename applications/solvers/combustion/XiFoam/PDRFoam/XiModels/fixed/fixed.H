/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::XiModels::fixed

Description
    Fixed value model for Xi. See \link XiModel.H \endlink for more details
    on flame wrinkling modelling.

SourceFiles
    fixed.C

\*---------------------------------------------------------------------------*/

#ifndef fixed_H
#define fixed_H

#include "XiModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                            Class fixed Declaration
\*---------------------------------------------------------------------------*/

class fixed
:
    public XiModel
{
public:

    //- Runtime type information
    TypeName("fixed");


    // Constructors

        //- Construct from components
        fixed
        (
            const dictionary& XiProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        fixed(const fixed&);


    //- Destructor
    virtual ~fixed();


    // Member Functions

        //- Correct the flame-wrinkling Xi
        virtual void correct()
        {}

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties);

        //- Write fields of the XiEq model
        virtual void writeFields()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fixed&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
