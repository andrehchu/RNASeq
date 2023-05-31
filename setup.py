from setuptools import setup, find_packages

setup(
    name='pySeqDiff',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'statsmodels',
        'scikit-learn',
        'argparse',
    ],
    entry_points='''
        [console_scripts]
        pySeqDiff=pySeqDiff.pySeqDiff:main
    ''',
)
