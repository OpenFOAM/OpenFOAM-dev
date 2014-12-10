{
    bool found = false;
    scalar random = rndGen.sample01<scalar>();
    while (!found && (n<99))
    {
        if (rrd_[n]>random)
        {
            found = true;
        }
        n++;
    }
    rNew = 0.04*n*rs;
}


