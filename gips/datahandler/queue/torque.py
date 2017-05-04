"""Torque-based backend for the GIPS work scheduler."""

# originally copied from:
# https://raw.githubusercontent.com/jdherman/hpc-submission-scripts/master/submit-pbs-loop.py

import os
from subprocess import check_output, Popen, PIPE, CalledProcessError

from gips import utils
from gips.datahandler.logger import Logger
from .. import queue

pbs_directives = [
    # for meaning of directives see
    #   http://docs.adaptivecomputing.com/torque/4-0-2/Content/topics/commands/qsub.htm
    '-N bsbreeder_job',
    '-A bsbreeder_test_contract:bsbreeder_test_task',
    '-k oe',
    '-j oe',
    '-V',
    '-l walltime=8:00:00',
]

import_block = """
import os
import datetime
import gippy
from gips.inventory import dbinv, orm
from gips.datahandler import worker
"""

setup_block = """
os.environ['GIPS_ORM'] = 'true'
gippy.Options.SetVerbose(4) # substantial verbosity for testing purposes
orm.setup()
"""

def get_job_name ():
    """Returns the current job ID, if one exists."""
    if 'PBS_JOBID' not in os.environ:
        raise queue.NoCurrentJobError('PBS_JOBID not set; no current job found')
    return os.environ['PBS_JOBID']


def is_job_alive (sched_id):
    try:
        result = check_output(['qstat', sched_id])
    except CalledProcessError as e:
        return False

    lines = result.split('\n')
    status = lines[2].split()[4]
    print status
    if status == 'C':
        return False
    else:
        return True;


def generate_script(operation, args_batch):
    """Produce torque script from given worker function and arguments."""
    lines = []
    lines.append('#!' + utils.settings().REMOTE_PYTHON) # shebang
    lines += ['#PBS ' + d for d in pbs_directives]      # #PBS directives
    try:
        outfile = utils.settings().TORQUE_OUTPUT
        lines.append('#PBS -o ' + outfile)
    except:
        pass
    lines.append(import_block)  # python imports, std lib, 3rd party, and gips
    lines.append(setup_block)   # config & setup code

    lines.append("print 'starting on `{}` job, sample arguments:'".format(operation))
    # TODO: this is risky. if args_batch has quotes in it, generates syntax error
    #       i switched to double quotes because repr seems to generate singles in most cases
    lines.append('print "{}"'.format(repr(args_batch[0])))

    # star of the show, the actual operation
    for args in args_batch:
        lines.append("worker.{}{}".format(operation, tuple(args)))

    return '\n'.join(lines) # stitch into single string & return


def submit(operation, args_ioi, batch_size=None, nproc=1, chain=False):
    """Submit jobs to the configured Torque system.

    Return value is a list of tuples, one for each batch:
        (job identifier, that batch's part of args_ioi, exit status of qsub, qsub's stdout, qsub's stderr)

    operation:  Defines which function will be performed, and must be one of
        'fetch', 'process', 'export', or 'postprocess'.
    args_ioi:  An iterable of iterables; each inner iterable gives the
        arguments to one call to the chosen function.
        The first item in the inner iterable needs to be the pk of the
        item to be fetched/queried/exported/etc.
    batch_size:  The work is divided among torque jobs; each job receives
        batch_size function calls to perform in a loop.  Leave None for one job
        that works the whole batch.
    nproc: number of processors to request
    chain: if True, chain batches to run in sequence
    """
    if batch_size is None:
        chunks = [args_ioi]
    else:
        chunks = list(utils.grouper(args_ioi, batch_size))

    # clean last chunk:  needed due to izip_longest padding chunks with Nones:
    chunks.append([i for i in chunks.pop() if i is not None])

    outcomes = []
    last_id = None  # used to chain jobs, if requested

    qsub_cmd = ['qsub']
    try:
        queue = '-q' + utils.settings().TORQUE_QUEUE
        qsub_cmd.append(queue)
    except:
        pass

    try:
        node = utils.settings().TORQUE_NODE
    except:
        node = 1
    qsub_cmd.append('-lnodes={}:ppn={}'.format(node, nproc))
    
    for chunk in chunks:
        job_script = generate_script(operation, chunk)

        if chain and last_id is not None:
            qsub = list(qsub_cmd)
            qsub.append('-Wdepend=afterany:{}'.format(last_id))
        else:
            qsub = qsub_cmd

        # open a pipe to the qsub command, then end job_string to it
        proc = Popen(qsub, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
        proc.stdin.write(job_script)
        out, err = proc.communicate()
        if proc.returncode != 0:
            utils.verbose_out("qsub failed:\n" + err, 1)
            sched_id = None
            ids = None
        else:
            sched_id = out.strip()
            if chain:
                last_id = sched_id
        outcomes.append((sched_id, chunk, proc.returncode, out, err))

        # if debug mode, print the generated torque script to a file
        if utils.get_setting('DH_DEBUG'):
            dfn = '{}-{}-job-script.py'.format(operation, sched_id)
            utils.verbose_out('dumping to ' + dfn, 1)
            with open(dfn, 'w') as dump_file:
                dump_file.write(job_script)

    # TODO confirm qsub exited 0 (raise otherwise)
    # TODO return best thing for checking on status
    # TODO log err someplace

    return outcomes
