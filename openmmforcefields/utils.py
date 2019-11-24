from pkg_resources import resource_filename
import contextlib, time

def get_ffxml_path():
    """
    Return the path where OpenMM ffxml forcefield files are stored in this package.

    Returns
    -------
    path : str
        The absolute path where OpenMM ffxml forcefield files are stored in this package
    """
    filename = resource_filename('openmmforcefields', 'ffxml')
    return filename

def get_data_filename(relative_path):
    """get the full path to one of the reference files shipped for testing

    in the source distribution, these files are in ``perses/data/*/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    relative_path : str
        name of the file to load (with respect to the openmoltools folder).

    Returns
    -------
    Absolute path to file

    """
    from pkg_resources import resource_filename
    import os

    fn = resource_filename('openmmforcefields', 'data/' + relative_path)

    if not os.path.exists(fn):
        raise ValueError("sorry! %s does not exist. if you just added it, you'll have to re-install" % fn)

    return fn

# =============================================================================
# BENCHMARKING UTILITIES
# =============================================================================

@contextlib.contextmanager
def time_it(task_name):
    """Context manager to log execution time of a block of code.

    Parameters
    ----------
    task_name : str
        The name of the task that will be reported.

    """
    timer = Timer()
    timer.start(task_name)
    yield timer  # Resume program
    timer.stop(task_name)
    timer.report_timing()


def with_timer(task_name):
    """Decorator that logs the execution time of a function.

    Parameters
    ----------
    task_name : str
        The name of the task that will be reported.

    """
    def _with_timer(func):
        @functools.wraps(func)
        def _wrapper(*args, **kwargs):
            with time_it(task_name):
                return func(*args, **kwargs)
        return _wrapper
    return _with_timer


class Timer(object):
    """A class with stopwatch-style timing functions.

    Examples
    --------
    >>> timer = Timer()
    >>> timer.start('my benchmark')
    >>> for i in range(10):
    ...     pass
    >>> elapsed_time = timer.stop('my benchmark')
    >>> timer.start('second benchmark')
    >>> for i in range(10):
    ...     for j in range(10):
    ...         pass
    >>> elsapsed_time = timer.stop('second benchmark')
    >>> timer.report_timing()

    """

    def __init__(self):
        self.reset_timing_statistics()

    def reset_timing_statistics(self, benchmark_id=None):
        """Reset the timing statistics.

        Parameters
        ----------
        benchmark_id : str, optional
            If specified, only the timings associated to this benchmark
            id will be reset, otherwise all timing information are.

        """
        if benchmark_id is None:
            self._t0 = {}
            self._t1 = {}
            self._completed = {}
        else:
            self._t0.pop(benchmark_id, None)
            self._t1.pop(benchmark_id, None)
            self._completed.pop(benchmark_id, None)

    def start(self, benchmark_id):
        """Start a timer with given benchmark_id."""
        self._t0[benchmark_id] = time.time()

    def stop(self, benchmark_id):
        try:
            t0 = self._t0[benchmark_id]
        except KeyError:
            logger.warning("Can't stop timing for {}".format(benchmark_id))
        else:
            self._t1[benchmark_id] = time.time()
            elapsed_time = self._t1[benchmark_id] - t0
            self._completed[benchmark_id] = elapsed_time
            return elapsed_time

    def partial(self, benchmark_id):
        """Return the elapsed time of the given benchmark so far."""
        try:
            t0 = self._t0[benchmark_id]
        except KeyError:
            logger.warning("Couldn't return partial timing for {}".format(benchmark_id))
        else:
            return time.time() - t0

    def report_timing(self, clear=True):
        """Log all the timings at the debug level.

        Parameters
        ----------
        clear : bool
            If True, the stored timings are deleted after being reported.

        Returns
        -------
        elapsed_times : dict
            The dictionary benchmark_id : elapsed time for all benchmarks.

        """
        for benchmark_id, elapsed_time in self._completed.items():
            logger.debug('{} took {:8.3f}s'.format(benchmark_id, elapsed_time))

        if clear is True:
            self.reset_timing_statistics()
