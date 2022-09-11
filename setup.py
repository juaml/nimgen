# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
# License: AGPL
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

DOWNLOAD_URL = 'https://github.com/juaml/nimgen/'
URL = 'https://juaml.github.io/nimgen'

setuptools.setup(
    name='nimgen',
    author='Applied Machine Learning',
    author_email='sami.hamdan@fz-juelich.de',
    description='FZJ AML NIMGEN Library ',
    long_description=long_description,
    entry_points={
        "console_scripts": ["nimgen = nimgen.nimgen:main"]
    },
    long_description_content_type='text/markdown',
    url=URL,
    download_url=DOWNLOAD_URL,
    packages=setuptools.find_packages(),
    zip_safe=False,
    classifiers=['Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'License :: OSI Approved',
                 'Programming Language :: Python',
                 'Topic :: Software Development',
                 'Topic :: Scientific/Engineering',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS',
                 'Programming Language :: Python :: 3'],
    project_urls={
        'Documentation': URL,
        'Source': DOWNLOAD_URL,
        'Tracker': f'{DOWNLOAD_URL}issues/',
    },
    install_requires=['numpy>=1.19.1',
                      'pandas>=1.1.2',
                      'abagen'],
    python_requires='>=3.6',
    package_data={'': ['*.r', '*.R']},
    include_package_data=True,
)
