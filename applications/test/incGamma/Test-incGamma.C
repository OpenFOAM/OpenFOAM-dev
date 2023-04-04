/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

Application
    Test-incGamma

Description
    This checks the consistency of incGamma and invIncGamma. It writes files
    with columns containing values for different 'a' parameters. The
    incGammaRatio_P.xy and invIncGammaRatio_P.xy files should generate the same
    lines when the plot axis order is reversed. E.g., in gnuplot:

        plot "incGammaRatio_P.xy" us 1:5 w l, \
             "invIncGammaRatio_P.xy" us 5:1 w p

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "mathematicalConstants.H"
#include "rawSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace Foam::constant::mathematical;

scalar grade(const scalar f)
{
    return (Foam::cos(pi*(1 - f)) + 1)/2;
}

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    const scalarField as({0.25, 0.5, 1, 2, 4});

    static const scalar nXs = 10000;
    static const scalar xMax = 10;
    static const scalar nPs = 100;

    wordList aNames(as.size());
    forAll(as, ai)
    {
        aNames[ai] = name(as[ai]);
    }

    scalarField xs(nXs);
    forAll(xs, xi)
    {
        const scalar f = grade(grade(scalar(xi)/(nXs - 1)));
        xs[xi] = (1 - f)*small + f*xMax;
    }

    scalarField Ps(nPs);
    forAll(Ps, Pi)
    {
        const scalar f = grade(grade(scalar(Pi)/(nPs - 1)));
        Ps[Pi] = (1 - f)*rootVSmall + f*(1 - 1e-4);
    }

    #define DeclareTypePXs(Type, nullArg) \
        PtrList<Field<Type>> Type##Ps(as.size()); \
        PtrList<Field<Type>> Type##PErrors(as.size()); \
        PtrList<Field<Type>> Type##Xs(as.size());
    FOR_ALL_FIELD_TYPES(DeclareTypePXs);
    #undef DeclareTypePXs

    forAll(as, ai)
    {
        scalarPs.set(ai, new scalarField(nXs));
        scalarPErrors.set(ai, new scalarField(nXs));
        forAll(xs, xi)
        {
            const scalar P = incGammaRatio_P(as[ai], xs[xi]);
            const scalar PStar =
                incGammaRatio_P
                (
                    as[ai],
                    invIncGammaRatio_P
                    (
                        as[ai],
                        min(max(P, small), 1 - small)
                    )
                );
            scalarPs[ai][xi] = P;
            scalarPErrors[ai][xi] = mag(P - PStar);
        }
    }
    forAll(as, ai)
    {
        scalarXs.set(ai, new scalarField(nPs));
        forAll(Ps, Pi)
        {
            scalarXs[ai][Pi] = invIncGammaRatio_P(as[ai], Ps[Pi]);
        }
    }

    rawSetWriter(IOstream::ASCII, IOstream::UNCOMPRESSED).write
    (
        ".",
        "incGammaRatio_P",
        coordSet(true, "x", xs),
        aNames
        #define TypePsParameter(Type, nullArg) , Type##Ps
        FOR_ALL_FIELD_TYPES(TypePsParameter)
        #undef TypePsParameter
    );
    rawSetWriter(IOstream::ASCII, IOstream::UNCOMPRESSED).write
    (
        ".",
        "incGammaRatio_PError",
        coordSet(true, "x", xs),
        aNames
        #define TypePErrorsParameter(Type, nullArg) , Type##PErrors
        FOR_ALL_FIELD_TYPES(TypePErrorsParameter)
        #undef TypePErrorsParameter
    );
    rawSetWriter(IOstream::ASCII, IOstream::UNCOMPRESSED).write
    (
        ".",
        "invIncGammaRatio_P",
        coordSet(true, "P", Ps),
        aNames
        #define TypeXsParameter(Type, nullArg) , Type##Xs
        FOR_ALL_FIELD_TYPES(TypeXsParameter)
        #undef TypeXsParameter
    );

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
