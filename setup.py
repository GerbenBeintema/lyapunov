from setuptools import setup, find_namespace_packages

with open('requirements.txt') as f:
    install_requires = [line for line in f]

name = 'lyapunov'
packages = [a for a in find_namespace_packages(where='.') if a[:len(name)]==name]
setup(name = name,
      version = '0.1.0',
      description = 'lyapunov: An easy-going maximum lyapunov exponent computer for continuous systems',
      author = 'Gerben Beintema',
      author_email = 'g.i.beintema@tue.nl',
      license = 'BSD 3-Clause License',
      python_requires = '>=3.6',
      packages=packages,
      install_requires = install_requires,
    )