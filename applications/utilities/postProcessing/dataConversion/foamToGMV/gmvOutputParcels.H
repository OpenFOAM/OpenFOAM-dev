gmvFile << "tracers " << particles.size() << nl;
forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().position().x() << " ";
}
gmvFile << nl;

forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().position().y() << " ";
}
gmvFile << nl;

forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().position().z() << " ";
}
gmvFile << nl;

gmvFile << "U" << nl;
forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().velocity().x() << " ";
}
gmvFile << nl;

gmvFile << "V" << nl;
forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().velocity().y() << " ";
}
gmvFile << nl;

gmvFile << "W" << nl;
forAllConstIter(discretePhase, particles, iter)
{
{
    gmvFile << iter().velocity().z() << " ";
}
gmvFile << nl;

gmvFile << "Diam" << nl;
forAllConstIter(discretePhase, particles, iter)
{
    gmvFile << iter().d() << " ";
}

gmvFile << "endtrace"<< nl;
