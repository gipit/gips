"""Interface for task queues; currently supporting RQ and TORQUE.

The API has three functions, each of which is a pass-through to specific
code that talks to RQ or TORQUE, depending on the GIPS setting
TASK_QUEUE.
"""

from gips import utils


valid_tq_settings = ('torque', 'rq')

_tq_setting = None


class NoCurrentJobError(RuntimeError):
    """When code expects to be "in" a running job, but is not."""
    pass


def tq_setting():
    """Return the value of gips setting TASK_QUEUE (a string).

    Uses a global to avoid re-reading settings multiple times.  Add to
    queue.valid_tq_settings to add support for a new queue type."""
    global _tq_setting
    if _tq_setting is None:
        _tq_setting = utils.get_setting('TASK_QUEUE', default='torque')
        if _tq_setting not in valid_tq_settings:
            raise ValueError('No such TASK_QUEUE `{}`'.format(value))
    return _tq_setting


def get_queue_module():
    """Return the appropriate module based on GIPS setting TASK_QUEUE."""
    if tq_setting() == 'torque':
        from . import torque
        return torque
    if tq_setting() == 'rq':
        from . import rq
        return rq

    raise RuntimeError("Unreachable line reached; valid_tq_settings need updating?")


def get_job_name(*args, **kwargs):
    """Return the job's name as given by its queueing system."""
    # need to return an ascii string because, way down the line, gippy can't handle unicode
    return get_queue_module().get_job_name(*args, **kwargs).encode('ascii')

def is_job_alive(*args, **kwargs):
    """Returns True only if the job exists and is not yet completed.

    This means that a 'living' job can running, or can be in a work queue but
    not yet started.
    """
    return get_queue_module().is_job_alive(*args, **kwargs)


def submit(operation, call_signatures, batch_size=None, nproc=1, chain=False):
    """Not all of these arguments are supported by all task modules.

    See the task modules for details.
    """
    valid_ops = ('query', 'fetch', 'process', 'export', 'export_and_aggregate')
    assert operation in valid_ops, '{} is not a valid operation'.format(operation)

    qm = get_queue_module()

    # here down is full of hacks to avoid fixing hacks in torque and scheduler
    if tq_setting() == 'torque':
        outcomes = qm.submit(operation, call_signatures, batch_size, nproc, chain)
        return outcomes

    if tq_setting() == 'rq':
        task_ids = qm.submit(operation, call_signatures, chain)
        # scheduler expects batches of signatures per task ID, hence [cs]
        digestible_by_scheduler = zip(task_ids, [[cs] for cs in call_signatures])
        return digestible_by_scheduler

    raise RuntimeError("Unreachable line reached; valid_tq_settings need updating?")
