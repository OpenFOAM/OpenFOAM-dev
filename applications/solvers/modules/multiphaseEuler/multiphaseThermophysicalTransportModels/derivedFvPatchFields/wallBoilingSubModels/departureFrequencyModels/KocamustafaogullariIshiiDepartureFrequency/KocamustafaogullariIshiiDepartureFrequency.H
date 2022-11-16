/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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
    Foam::wallBoilingModels::departureFrequencyModels::
    KocamustafaogullariIshiiDepartureFrequency

Description
    Correlation for bubble departure frequency.

    Reference:
    \verbatim
        Kocamustafaogullari, G., & Ishii, M. (1995).
        Foundation of the interfacial area transport equation and its closure
        relations.
        International Journal of Heat and Mass Transfer, 38(3), 481-493.
    \endverbatim

SourceFiles
    KocamustafaogullariIshiiDepartureFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef KocamustafaogullariIshiiDepartureFrequency_H
#define KocamustafaogullariIshiiDepartureFrequency_H

#include "departureFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureFrequencyModels
{

/*---------------------------------------------------------------------------*\
         Class KocamustafaogullariIshiiDepartureFrequency Declaration
\*---------------------------------------------------------------------------*/

class KocamustafaogullariIshiiDepartureFrequency
:
    public departureFrequencyModel
{
    // Private data

        //- Model parameter
        scalar Cf_;


public:

    //- Runtime type information
    TypeName("KocamustafaogullariIshii");


    // Constructors

        //- Construct from a dictionary
        KocamustafaogullariIshiiDepartureFrequency(const dictionary& dict);

        //- Copy construct
        KocamustafaogullariIshiiDepartureFrequency
        (
            const KocamustafaogullariIshiiDepartureFrequency& model
        );

        //- Construct and return a clone
        virtual autoPtr<departureFrequencyModel> clone() const
        {
            return autoPtr<departureFrequencyModel>
            (
                new KocamustafaogullariIshiiDepartureFrequency(*this)
            );
        }


    //- Destructor
    virtual ~KocamustafaogullariIshiiDepartureFrequency();


    // Member Functions

        //- Calculate and return the bubble departure frequency
        virtual tmp<scalarField> fDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L,
            const scalarField& dDep
        ) const;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace departureFrequencyModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
