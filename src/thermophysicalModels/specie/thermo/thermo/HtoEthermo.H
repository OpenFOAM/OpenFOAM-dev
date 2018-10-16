inline scalar Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T) - EquationOfState::CpMCv(p, T);
}

inline scalar Es
(
    const scalar p,
    const scalar T
) const
{
    return Hs(p, T) - p/EquationOfState::rho(p, T);
}

inline scalar Ea
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T) - p/EquationOfState::rho(p, T);
}
