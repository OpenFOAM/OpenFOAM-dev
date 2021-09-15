/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::TroeFallOffFunction::TroeFallOffFunction
(
    const scalar alpha,
    const scalar Tsss,
    const scalar Ts,
    const scalar Tss
)
:
    alpha_(alpha),
    Tsss_(Tsss),
    Ts_(Ts),
    Tss_(Tss)
{}


inline Foam::TroeFallOffFunction::TroeFallOffFunction(const dictionary& dict)
:
    alpha_(dict.lookup<scalar>("alpha")),
    Tsss_(dict.lookup<scalar>("Tsss")),
    Ts_(dict.lookup<scalar>("Ts")),
    Tss_(dict.lookup<scalar>("Tss"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::TroeFallOffFunction::operator()
(
    const scalar T,
    const scalar Pr
) const
{
    const scalar logPr = log10(max(Pr, small));

    const scalar Fcent =
        (1 - alpha_)*exp(-T/Tsss_)
      + alpha_*exp(-T/Ts_)
      + exp(-Tss_/T);

    const scalar logFcent = log10(max(Fcent, small));

    const scalar c = -0.4 - 0.67*logFcent;
    static const scalar d = 0.14;
    const scalar n = 0.75 - 1.27*logFcent;

    const scalar x1 = n - d*(logPr + c);
    const scalar x2 = (logPr + c)/x1;
    const scalar x3 = 1 + sqr(x2);
    const scalar x4 = logFcent/x3;

    return pow(10, x4);
}


inline Foam::scalar Foam::TroeFallOffFunction::ddT
(
    const scalar T,
    const scalar Pr,
    const scalar F
) const
{
    static const scalar logTen = log(scalar(10));

    const scalar logPr = log10(max(Pr, small));

    const scalar Fcent =
        (1 - alpha_)*exp(-T/Tsss_)
      + alpha_*exp(-T/Ts_)
      + exp(-Tss_/T);
    const scalar dFcentdT =
      - (1 - alpha_)/Tsss_*exp(-T/Tsss_)
      - alpha_/Ts_*exp(-T/Ts_)
      + Tss_/sqr(T)*exp(-Tss_/T);

    const scalar logFcent = log10(max(Fcent, small));
    const scalar dlogFcentdT = Fcent >= small ? dFcentdT/Fcent/logTen : 0;

    const scalar c = -0.4 - 0.67*logFcent;
    const scalar dcdT = -0.67*dlogFcentdT;
    static const scalar d = 0.14;
    const scalar n = 0.75 - 1.27*logFcent;
    const scalar dndT = - 1.27*dlogFcentdT;

    const scalar x1 = n - d*(logPr + c);
    const scalar dx1dT = dndT - d*dcdT;
    const scalar x2 = (logPr + c)/x1;
    const scalar dx2dT = (dcdT - x2*dx1dT)/x1;
    const scalar x3 = 1 + sqr(x2);
    const scalar dx3dT = 2*x2*dx2dT;
    const scalar x4 = logFcent/x3;
    const scalar dx4dT = (dlogFcentdT - x4*dx3dT)/x3;

    return logTen*F*dx4dT;
}


inline Foam::scalar Foam::TroeFallOffFunction::ddPr
(
    const scalar T,
    const scalar Pr,
    const scalar F
) const
{
    static const scalar logTen = log(scalar(10));

    const scalar logPr = log10(max(Pr, small));
    const scalar dlogPrdPr = Pr >= small ? 1/(logTen*Pr) : 0;

    const scalar Fcent =
        (1 - alpha_)*exp(-T/Tsss_)
      + alpha_*exp(-T/Ts_)
      + exp(-Tss_/T);

    const scalar logFcent = log10(max(Fcent, small));

    const scalar c = -0.4 - 0.67*logFcent;
    static const scalar d = 0.14;
    const scalar n = 0.75 - 1.27*logFcent;

    const scalar x1 = n - d*(logPr + c);
    const scalar dx1dPr = -d*dlogPrdPr;
    const scalar x2 = (logPr + c)/x1;
    const scalar dx2dPr = (dlogPrdPr - x2*dx1dPr)/x1;
    const scalar x3 = 1 + sqr(x2);
    const scalar dx3dPr = 2*x2*dx2dPr;
    const scalar x4 = logFcent/x3;
    const scalar dx4dPr = -x4*dx3dPr/x3;

    return logTen*F*dx4dPr;
}


inline void Foam::TroeFallOffFunction::write(Ostream& os) const
{
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "Tsss", Tsss_);
    writeEntry(os, "Ts", Ts_);
    writeEntry(os, "Tss", Tss_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::TroeFallOffFunction& tfof
)
{
    tfof.write(os);
    return os;
}


// ************************************************************************* //
