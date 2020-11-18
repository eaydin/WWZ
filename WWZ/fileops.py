from WWZ.exceptions import TimeSeriesInconsistent


def read_timeseries(fp):
    """
    Read time series from file

    Parameters
    ----------
    fp: File pointer

    Returns
    -------
    tuple
        The first element is a list of time values. The second element is a list of magnitude values.

    """

    time_series = []
    magnitude = []

    for line in fp:
        # Check if it's a comment line
        if line.strip().startswith('%') or line.strip().startswith('#'):
            continue
        else:
            time_series.append(float(line.split()[0]))
            magnitude.append(float(line.split()[1]))
    fp.close()

    if len(time_series) != len(magnitude):
        raise TimeSeriesInconsistent("The number of elements in the timeseries don't match. Check your timeseries data.")

    return time_series, magnitude
