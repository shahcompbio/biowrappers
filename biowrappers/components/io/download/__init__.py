import pypeliner
from pypeliner.workflow import Workflow

import tasks


def create_download_workflow(url, file_name):

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('url'),
        value=url
    )

    workflow.transform(
        name='download',
        ctx={'local': True},
        func=tasks.download_from_url,
        args=(
            pypeliner.managed.TempInputObj('url'),
            pypeliner.managed.OutputFile(file_name)
        )
    )

    return workflow
