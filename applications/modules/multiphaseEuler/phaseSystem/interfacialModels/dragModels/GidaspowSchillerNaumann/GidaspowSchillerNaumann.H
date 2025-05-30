/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::dragModels::GidaspowSchillerNaumann

Description
    Gidaspow, Schiller and Naumann drag model

    References:
    \verbatim
        Enwald, H., Peirano, E., & Almstedt, A. E. (1996).
        Eulerian two-phase flow theory applied to fluidisation.
        International Journal of Multiphase Flow, 22, 21-66.

        This is identical to the Wen and Yu, Rowe model Table 3.6 p.56  in
        van Wachem, B. G. M. (2000).
        Derivation, implementation, and validation of computer simulation models
        for gas-solid fluidised beds.
        PhD Thesis, TU Delft.
    \endverbatim

SourceFiles
    GidaspowSchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowSchillerNaumann_H
#define GidaspowSchillerNaumann_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                   Class GidaspowSchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class GidaspowSchillerNaumann
:
    public dispersedDragModel
{

public:

    //- Runtime type information
    TypeName("GidaspowSchillerNaumann");


    // Constructors

        //- Construct from a dictionary and an interface
        GidaspowSchillerNaumann
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~GidaspowSchillerNaumann();


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
