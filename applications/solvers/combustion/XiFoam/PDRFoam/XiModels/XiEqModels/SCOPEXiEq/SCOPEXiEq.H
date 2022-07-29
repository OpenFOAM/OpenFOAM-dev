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
    Foam::XiEqModels::SCOPEXiEq

Description
    Simple SCOPEXiEq model for XiEq based on SCOPEXiEqs correlation
    with a linear correction function to give a plausible profile for XiEq.
    See \link SCOPELaminarFlameSpeed.H \endlink for details on the SCOPE laminar
    flame speed model.

SourceFiles
    SCOPEXiEq.C

\*---------------------------------------------------------------------------*/

#ifndef SCOPEXiEq_H
#define SCOPEXiEq_H

#include "XiEqModel.H"
#include "SCOPELaminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class SCOPEXiEq Declaration
\*---------------------------------------------------------------------------*/

class SCOPEXiEq
:
    public XiEqModel
{
    // Private Data

        // Model constant
        scalar XiEqCoef_;

        // Model constant
        scalar XiEqExp_;

        // Model constant
        scalar lCoef_;

        //- Minimum Su
        dimensionedScalar SuMin_;

        //- Schelkin effect Model constant
        scalar uPrimeCoef_;

        //- Use sub-grid Schelkin effect
        bool subGridSchelkin_;

        //- The SCOPE laminar flame speed model used to obtain the
        //  Marstein number.  Note: the laminar flame speed need not be
        //  obtained form the same model.
        laminarFlameSpeedModels::SCOPE MaModel;


public:

    //- Runtime type information
    TypeName("SCOPEXiEq");


    // Constructors

        //- Construct from components
        SCOPEXiEq
        (
            const dictionary& XiEqProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        SCOPEXiEq(const SCOPEXiEq&);


    //- Destructor
    virtual ~SCOPEXiEq();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SCOPEXiEq&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
