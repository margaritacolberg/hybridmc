import os
from setuptools import setup, find_packages

# Define the directory where your executables are located
scripts_dir = 'hybridmc/py_bin'

# List all the Python scripts in that directory except __init__.py which is empty
scripts = [os.path.join(scripts_dir, f) for f in os.listdir(scripts_dir) if f != '__init__.py' and f.endswith('.py')]

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="hybridmc",
    version=os.environ.get('GIT_DESCRIBE_TAG', '0.1.0').strip('v'),
    author="Vigneshwar Rajesh",
    author_email="vignesh.rajesh@mail.utoronto.ca",
    description="A project using c++ for Monte Carlo sampling based simulations.",
    long_description="",
    zip_safe=False,
    python_requires=">=3.11",
    packages=find_packages(),
    scripts=scripts
)
