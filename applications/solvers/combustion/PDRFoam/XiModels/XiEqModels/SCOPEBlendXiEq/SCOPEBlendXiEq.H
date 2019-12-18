/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::XiEqModels::SCOPEBlend

Description
    Simple SCOPEBlendXiEq model for XiEq based on SCOPEXiEqs correlation
    with a linear correction function to give a plausible profile for XiEq.
    See @link SCOPELaminarFlameSpeed.H @endlink for details on the SCOPE
    laminar flame speed model.

SourceFiles
    SCOPEBlend.C

\*---------------------------------------------------------------------------*/

#ifndef SCOPEBlendXiEq_H
#define SCOPEBlendXiEq_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class SCOPEBlend Declaration
\*---------------------------------------------------------------------------*/

class SCOPEBlend
:
    public XiEqModel
{
    // Private Data

        //- Low turbulence intensity equilibrium Xi model
        autoPtr<XiEqModel> XiEqModelL_;

        //- High turbulence intensity equilibrium Xi model
        autoPtr<XiEqModel> XiEqModelH_;


public:

    //- Runtime type information
    TypeName("SCOPEBlend");


    // Constructors

        //- Construct from components
        SCOPEBlend
        (
            const dictionary& XiEqProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        SCOPEBlend(const SCOPEBlend&);


    //- Destructor
    virtual ~SCOPEBlend();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties)
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SCOPEBlend&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
