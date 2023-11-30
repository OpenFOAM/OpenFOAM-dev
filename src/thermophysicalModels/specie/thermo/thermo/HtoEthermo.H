inline scalar Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T) - EquationOfState::CpMCv(p, T);
}

inline scalar es
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) - p/EquationOfState::rho(p, T);
}

inline scalar ea
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - p/EquationOfState::rho(p, T);
}
