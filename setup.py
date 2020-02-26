from setuptools import setup, find_packages

setup(
    name="netcdf-swan",
    version="0.0.1",
    author="Samuel Johnson",
    author_email="samdljohnson@gmail.com",
    url="https://github.com/meracan/netcdf-swan",
    packages=find_packages('src', exclude=["tests", "tests.*"]),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'scipy',
        'numpy',
        #'netcdf4',
    ],
    python_requires='>=3.6',
)
