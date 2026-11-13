"""Set the GitHub repository description, website, and topics after gh auth login."""
import os
import subprocess

REPOSITORY = 'alalbiol/iberian-pines-ifn-validation'
DESCRIPTION = ('Reproducible biomass forecasts for five Iberian pines, '
               'validated against Spanish National Forest Inventories (IFN3/IFN4).')
HOMEPAGE = 'https://www.eotdl.com/datasets/iberian-pines-ifn-validation'
TOPICS = ['forestry', 'biomass', 'forest-inventory', 'pine', 'ifn3', 'ifn4',
          'reproducible-research', 'python']


def main():
    env = {k: v for k, v in os.environ.items() if k not in ('GH_TOKEN', 'GITHUB_TOKEN', 'GH_DEBUG')}
    subprocess.run(['gh', 'auth', 'status', '--hostname', 'github.com'], env=env, check=True)
    command = ['gh', 'repo', 'edit', REPOSITORY, '--description', DESCRIPTION, '--homepage', HOMEPAGE]
    for topic in TOPICS:
        command.extend(['--add-topic', topic])
    subprocess.run(command, env=env, check=True)
    print(f'Updated About fields for {REPOSITORY}.')


if __name__ == '__main__':
    main()
