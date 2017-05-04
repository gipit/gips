import pytest

from gips.inventory.dbinv import models

from gips.datahandler import worker

@pytest.mark.django_db
def t_worker_export_and_aggregate(mocker):
    # mocks
    m_export    = mocker.patch.object(worker, '_export')
    m_aggregate = mocker.patch.object(worker, '_aggregate')
    m_makedirs  = mocker.patch.object(worker.os, 'makedirs')
    m_rmtree    = mocker.patch.object(worker.shutil, 'rmtree')

    catalog_entry = models.DataVariable.objects.create(product='whatever')
    job = models.Job.objects.create(variable=catalog_entry, status='pp-scheduled',
                                    spatial='0', temporal='0')
    # call
    worker.export_and_aggregate(job.pk, {})

    # assertions
    job.refresh_from_db()
    assert job.status == 'complete'
    for m in m_export, m_aggregate, m_makedirs, m_rmtree:
        m.assert_called_once() # TODO check arguments
