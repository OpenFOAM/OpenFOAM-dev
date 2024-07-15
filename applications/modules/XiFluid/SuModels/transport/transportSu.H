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
    Foam::SuModels::transport

Description
    Transport model for the strained laminar flame speed

    References:
    \verbatim
        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

SourceFiles
    transport.C

\*---------------------------------------------------------------------------*/

#ifndef transportSu_H
#define transportSu_H

#include "unstrained.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SuModels
{

/*---------------------------------------------------------------------------*\
                          Class transport Declaration
\*---------------------------------------------------------------------------*/

class transport
:
    public unstrained
{
    // Private Data

        //- Laminar flame extinction strain-rate
        dimensionedScalar sigmaExt_;

        //- Minimum laminar flame-speed allowed for numerical stability
        dimensionedScalar SuMin_;

        //- Maximum laminar flame-speed allowed for numerical stability
        dimensionedScalar SuMax_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("transport");


    // Constructors

        //- Construct from components
        transport
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        transport(const transport&) = delete;


    //- Destructor
    virtual ~transport();


    // Member Functions

        //- Correct the laminar flame speed
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SuModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
